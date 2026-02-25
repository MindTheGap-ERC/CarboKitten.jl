# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
import ...Production: production_profile, is_benthic, is_pelagic, capped_production, NoProduction

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production

# ~/~ begin <<docs/src/components/production.md#production-input>>[init]
@kwdef struct Input <: AbstractInput
    insolation
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
        elseif is_benthic(p)
            set_attribute(output, "facies$(i)/type", "benthic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        end
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]

    p(state::AbstractState, wd::AbstractMatrix) =
        capped_production.(production_rates[:, na, na], insolation_func(state), wd[na, :, :], dt)

    p(state::AbstractState) = p(state, w(state))

    return p
end

end
# ~/~ end
