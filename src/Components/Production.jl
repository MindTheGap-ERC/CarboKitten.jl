# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using HDF5
using Logging

export production_rate, uniform_production

@kwdef struct Facies <: AbstractFacies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
end

@kwdef struct Input <: AbstractInput
    insolation
end

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

function write_header(fid, input::AbstractInput)
    if input.insolation isa Quantity
        fid["input"]["insolation"] = fill(input.insolation |> in_units_of(u"W/m^2"), input.time.steps)
    elseif input.insolation isa AbstractVector
        fid["input"]["insolation"] = input.insolation |> in_units_of(u"W/m^2")
    else
        t = write_times(input)[1:end-1]
        fid["input"]["insolation"] = input.insolation.(t) |> in_units_of(u"W/m^2")
    end

    # fid["input"]["insolation"] = insolation(input) |> in_units_of(u"W/m^2")
    for (i, f) in enumerate(input.facies)
        attr = attributes(fid["input/facies$(i)"])
        attr["maximum_growth_rate"] = f.maximum_growth_rate |> in_units_of(u"m/Myr")
        attr["extinction_coefficient"] = f.extinction_coefficient |> in_units_of(u"m^-1")
        attr["saturation_intensity"] = f.saturation_intensity |> in_units_of(u"W/m^2")
    end
end

# ~/~ begin <<docs/src/components/production.md#component-production-rate>>[init]
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end
# ~/~ end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies

    return function (state::AbstractState)
        return production_rate.(
            insolation_func(state),
            facies[:, na, na],
            w(state)[na, :, :])
    end
end
end
# ~/~ end
