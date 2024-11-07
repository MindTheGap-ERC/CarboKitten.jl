# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using HDF5

export production_rate, uniform_production

@kwdef struct Facies <: AbstractFacies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
end

@kwdef struct Input <: AbstractInput
    insolation::Intensity
end

function write_header(fid, input::AbstractInput)
    attributes(fid["input"])["insolation"] = input.insolation |> in_units_of(u"W/m^2")
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
    return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end
# ~/~ end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]

    return function (state::AbstractState)
        return production_rate.(
            input.insolation,
            input.facies[:, na, na],
            w(state)[na, :, :])
    end
end
end
# ~/~ end
