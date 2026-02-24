# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using HDF5
using QuadGK
using Interpolations
using Logging

export production_rate, capped_production, uniform_production, benthic_production, pelagic_production

# ~/~ begin <<docs/src/components/production.md#component-production-rate>>[init]
function benthic_production(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end

production_rate(i, f, w) = benthic_production(i, f, w)

function capped_production(f, insolation, water_depth, dt)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(insolation, water_depth))
    return min(max(0.0u"m", water_depth), p * dt)
end
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#pelagic-production>>[init]
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#production-profile>>[init]
"""
    production_profile(input::AbstractInput, p)

Given an input production configuration, returns a function of
insolation (in W/m^2) and water depth (in m) to production (in m/Myr).

The default implementation is the identity function; it assumes
that its parameter is a callable object with the correct behaviour.

## Example
If you want to implement your own production profile, define a new type,

    struct MyProduction <: AbstractProduction
        ...
    end

and implement this method for that type,

    import CarboKitten.Components.Production: production_profile

    production_profile(input::AbstractInput, p::MyProduction) =
        function (insolation, water_depth)
            ...
        end
"""
production_profile(input::AbstractInput, p) = p

"""
    is_benthic(obj)

Predicate to determine if a facies or production spec is benthic.
Defaults to `false`.
"""
function is_benthic end

"""
    is_pelagic(obj)

Predicate to determine if a facies or production spec is pelagic.
Defaults to `false`.
"""
function is_pelagic end

struct NoProduction <: AbstractProduction
end

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::Rate = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::Intensity = 1.0u"W/m^2"
end

production_profile(input::AbstractInput, p::BenthicProduction) = (I, w) -> benthic_production(I, p, w)
is_benthic(::BenthicProduction) = true
is_pelagic(::BenthicProduction) = false

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr^-1") = 0.0u"Myr^-1"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::Intensity = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

production_profile(input::AbstractInput, p::PelagicProduction) = (I, w) -> pelagic_production_lookup(I, p, w)
is_benthic(::PelagicProduction) = false
is_pelagic(::PelagicProduction) = true

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#production-input>>[init]
@kwdef struct Input <: AbstractInput
    insolation
end
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
# ~/~ begin <<docs/src/components/production.md#production-lookup>>[init]
function pelagic_production_lookup(input::AbstractInput, prod::PelagicProduction)
    insolation_param = input.insolation
    depth_grid = LinRange(0.0, prod.maximum_production_depth, prod.table_size[2])

    if insolation_param isa Quantity
        production_values = [pelagic_production(insolation_param, prod, w) for w in depth_grid]
        interpolated = linear_interpolation(depth_grid, production_values, extrapolation_bc = Line())
        return (I0, w) -> interpolated(w)
    end

    insolation_vec = if insolation_param isa AbstractVector
        insolation_param
    else
        t = time_axis(input)
        insolation_param.(t)
    end
    I_min, I_max = extrema(insolation_vec)

    insolation_grid = LinRange(I_min, I_max, prod.table_size[1])
    production_values = [pelagic_production(I, prod, w) for I in insolation_grid, w in depth_grid]
    return linear_interpolation((insolation_grid, depth_grid), production_values, extrapolation_bc = Line())
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
