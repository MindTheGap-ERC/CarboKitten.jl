# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using ...Production: NoProduction, InterpolatedProduction,
    AbstractProductionModifier, MultiplyProduction, production_factor
import ...Production: production_profile, is_benthic, is_pelagic, capped_production

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production
export AbstractProductionModifier, MultiplyProduction, InterpolatedProduction

# ~/~ begin <<docs/src/components/production.md#production-input>>[init]
@kwdef struct Input <: AbstractInput
    insolation
    production_modifiers = []
end

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#production-insolation>>[init]
function insolation(input::AbstractInput)
    insolation = input.insolation
    tprop = input.time
    if insolation isa Quantity
        return s -> insolation
    end
    if insolation isa AbstractVector
        @info "Reading insolation from a table"
        return s -> insolation[s.step+1]
    end
    @info "Reading insolation from a function"
    function (s::AbstractState)
        t = time(tprop, s)
        return insolation(t)
    end
end
# ~/~ end

function write_header(input::AbstractInput, output::AbstractOutput)
    if input.insolation isa Quantity
        set_attribute(output, "insolation",
            fill(input.insolation |> in_units_of(u"W/m^2"), input.time.steps))
    elseif input.insolation isa AbstractVector
        set_attribute(output, "insolation",
            input.insolation |> in_units_of(u"W/m^2"))
    else
        t = write_times(input)[1:end-1]
        set_attribute(output, "insolation",
            input.insolation.(t) |> in_units_of(u"W/m^2"))
    end

    for (i, f) in enumerate(input.facies)
        p = f.production
        if is_pelagic(p)
            set_attribute(output, "facies$(i)/type", "pelagic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"1/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        elseif p isa InterpolatedProduction
            set_attribute(output, "facies$(i)/type", "interpolated")
            set_attribute(output, "facies$(i)/maximum_production", p.maximum_production |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/depth_knots", p.depth_knots .|> in_units_of(u"m"))
            set_attribute(output, "facies$(i)/multipliers", collect(Float64, p.multipliers))
        elseif is_benthic(p)
            set_attribute(output, "facies$(i)/type", "benthic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        end
    end

    for (idx, m) in enumerate(input.production_modifiers)
        prefix = "production_modifiers/m$(idx)"
        set_attribute(output, "$(prefix)/kind", "MultiplyProduction")
        set_attribute(output, "$(prefix)/factor", m.factor)
        set_attribute(output, "$(prefix)/t_range",
            m.t_range isa Colon ? [NaN, NaN] : [ustrip(u"Myr", m.t_range[1]), ustrip(u"Myr", m.t_range[2])])
        set_attribute(output, "$(prefix)/facies",
            m.facies isa Colon ? [-1] : m.facies isa Int ? [m.facies] : collect(Int, m.facies))
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]

    modifiers = hasfield(typeof(input), :production_modifiers) ?
        input.production_modifiers : AbstractProductionModifier[]

    if isempty(modifiers)
        p(state::AbstractState, wd::AbstractMatrix) =
            capped_production.(production_rates[:, na, na], insolation_func(state), wd[na, :, :], dt)
        p(state::AbstractState) = p(state, w(state))
        return p
    end

    n_f = length(facies)
    get_time = time(input)
    function p_mod(state::AbstractState, wd::AbstractMatrix)
        t = get_time(state)
        factors = Float64[production_factor(modifiers, t, i) for i in 1:n_f]
        return capped_production.(production_rates[:, na, na], insolation_func(state),
                                  wd[na, :, :], dt, factors[:, na, na])
    end
    p_mod(state::AbstractState) = p_mod(state, w(state))
    return p_mod

end

end
# ~/~ end
