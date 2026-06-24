# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using ...Production: NoProduction, InterpolatedProduction, MultiplyProduction
import ...Production: production_profile, is_benthic, is_pelagic, is_interpolated,
    capped_production, insolation_curve

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production
export MultiplyProduction, InterpolatedProduction, NoProduction

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

function write_header(input::AbstractInput, output::AbstractOutput)
    # Insolation time series
    I_of_t = insolation_curve(input)
    t_write = write_times(input)[1:end-1]
    set_attribute(output, "insolation",
        [I_of_t(t) |> in_units_of(u"W/m^2") for t in t_write])

    # Generic 2D production table — no type-specific branches needed
    depth_grid = LinRange(0.0u"m", 200.0u"m", 200)
    for (i, f) in enumerate(input.facies)
        prof = production_profile(input, f.production)
        table = [prof(t, d) |> in_units_of(u"m/Myr")
                 for d in depth_grid, t in t_write]
        set_attribute(output, "facies$(i)/production_table", table)
        set_attribute(output, "facies$(i)/production_depth_axis",
            collect(depth_grid) .|> in_units_of(u"m"))
        set_attribute(output, "facies$(i)/production_time_axis",
            t_write .|> in_units_of(u"Myr"))
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]
    get_time = time(input)

    p(state::AbstractState, wd::AbstractMatrix) = begin
        t = get_time(state)
        capped_production.(production_rates[:, na, na], t, wd[na, :, :], dt)
    end
    p(state::AbstractState) = p(state, w(state))
    return p
end

end
# ~/~ end
