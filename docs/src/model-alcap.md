# Model with CA, Production and Active Layer transport (ALCAPS)

The following **S**edimentation model includes the Burgess 2013 **C**ellular **A**utomaton, Bosscher & Schlager 1992 **P**roduction curves and an **A**ctive **L**ayer transport model, based on Paola 1992, henceforth ALCAPS.

``` {.julia file=src/Model/ALCAPS.jl}
module ALCAPS

using Unitful
using HDF5

using ...BoundaryTrait: Shelf
using ...Config: Box, TimeProperties, axes
using ...SedimentStack: push_sediment!, pop_sediment!
using ...Burgess2013.CA: step_ca
using ...Burgess2013.Production: production_rate
using ...Transport.ActiveLayer: pde_stencil, Rate, Amount
using ...Utility: in_units_of

<<alcaps>>

end
```

## Input
It is convenient to define the `Rate` and `Amount` types.

``` {.julia #alcaps}
const Myr = u"Myr"
const m = u"m"

const PERIOD = 200.0u"kyr"
const AMPLITUDE = 4.0u"m"

@kwdef struct Facies
    viability_range::Tuple{Int,Int}
    activation_range::Tuple{Int,Int}

    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    # though units are in m, this in not an amount.
    # TODO: figure out what this unit means
    # values should be picked rather large, say 10km.
    diffusion_coefficient::typeof(1.0m)
end

const MODEL1 = [
    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 500u"m/Myr",
           extinction_coefficient = 0.8u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 1000u"m"),

    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 400u"m/Myr",
           extinction_coefficient = 0.1u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 5000u"m"),

    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 100u"m/Myr",
           extinction_coefficient = 0.005u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 10000u"m")
]

@kwdef struct Input
    box::Box              = Box{Shelf}(grid_size=(100, 50), phys_scale=150m)
    time::TimeProperties  = TimeProperties(𝚫t=0.001Myr, steps=1000, write_interval=1)
    ca_interval::Int      = 10

    bedrock_elevation     = (x, y) -> -x / 300.0  # (m, m) -> m
    sea_level             = t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate::Rate = 50.0m/Myr

    facies::Vector{Facies}    = MODEL1
    disintegration_rate::Rate = 500.0m/Myr   # same as maximum production rate
    insolation::typeof(1.0u"W/m^2") = 400.0u"W/m^2"

    sediment_buffer_size::Int     = 50
    sediment_buffer_grain::Amount = 0.5m
end
```

## Logic
It doesn't really matter if we pass around amounts or rates here. The equations solve the same.

``` {.julia #alcaps}
## Emitting a ModelFrame every iteration allows for
## inspecting the output of an entire run
## For 100x50x3 x 3 x 10000 x 8B ≈ 5GB
struct ModelFrame
    disintegration::Array{Amount,3}    # facies, x, y
    production::Array{Amount,3}
    deposition::Array{Amount,3}
end

mutable struct State
    time::typeof(1.0u"Myr")

    ca::Array{Int}
    ca_priority::Vector{Int}

    sediment_height::Array{Amount,2}   # x, y
    # sediment_buffer stores fractions, so no units
    sediment_buffer::Array{Float64,4}  # facies, x, y, z
end

function initial_state(input)
    sediment_height = zeros(Float64, input.box.grid_size...) * u"m"
    n_facies = length(input.facies)
    sediment_buffer = zeros(Float64, n_facies, input.box.grid_size..., input.sediment_buffer_size)

    ca = rand(0:n_facies, input.box.grid_size...)
    state = State(0.0u"Myr", ca, 1:n_facies, sediment_height, sediment_buffer)

    for _ = 1:20
        step_ca(input.box, input.facies)(state)
    end

    return state
end

function disintegration(input)
    max_h = input.disintegration_rate * input.Δt

    return function(state)
        h = min.(max_h, state.sediment_height)
        state.sediment_height .-= h
        return pop_sediment!(state.sediment_buffer, h)
    end
end

function production(input)
    p(f, w) = production_rate(input.insolation, input.facies[f], w)
    x, y = axes(input.box)
    w0 = .- input.bedrock_elevation(x, y')

    return function(state)
        water_depth = w0 .+ input.sea_level(state.time) .-
            state.sediment_height .+ input.subsidence_rate * state.time
        return p.(state.ca, water_depth)
    end
end

function transportation(input)
    n_facies = length(input.facies)
    x, y = axes(input.box)
    μ0 = input.bedrock_elevation.(x, y')
    # We always return this array
    transported_output = Array{Amount, 3}(undef, n_facies, input.box.grid_size...)
    stencils = [
        let stc = pde_stencil(input.box, f.diffusion_coefficient)
            (μ, p) -> @views stc(tuple.(μ, p[i,:,:]), transported_output[i,:,:])
        end for (i, f) in enumerate(input.facies) ]

    return function(state, active_layer::Array{Amount, 3})
        μ = state.sediment_height .+ μ0
        for stc in stencils
            stc(μ, active_layer)
        end

        return transported_output
    end
end

function run_model(input)
    state = initial_state(input)
    step_ca! = step_ca(input.box, input.facies)
    disintegrate! = disintegration(input)
    produce = production(input)
    transport = transportation(input)

    Channel{ModelFrame}() do ch
        for i in 1:input.time.steps
            if mod(i, input.ca_interval) == 0
                step_ca!(state)
            end

            p = produce(state)
            d = disintegrate!(state)

            active_layer = p.production .+ d.disintegration
            sediment = transport(state, active_layer)

            put!(ch, ModelFrame(d, p, sediment))

            push_sediment!(state.sediment_buffer, sediment ./ input.sediment_buffer_grain |> NoUnits)
            state.sediment_height .+= sum(sediment; dims=1)
            state.time += input.time.𝚫t
        end
    end
end

function main(input::Input, output::String)
    x, y = axes(input.box)

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x) |> in_units_of(u"m")
        gid["y"] = collect(y) |> in_units_of(u"m")
        gid["bedrock_elevation"] = input.bedrock_elevation.(x, y') |> in_units_of(u"m")
        gid["t"] = collect(0:(input.time.steps-1) .* input.time.Δt) |> in_units_of(u"Myr")

        attr = attributes(gid)
        attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
        attr["write_interval"] = input.time.write_interval
        attr["time_steps"] = input.time.steps
        attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")
        attr["n_facies"] = length(input.facies)

        n_facies = length(input.facies)
        ds_prod = create_dataset(fid, "production", datatype(Float64),
            dataspace(input.box.grid_size..., n_facies, input.time.steps),
            chunk=(input.box.grid_size..., n_facies, 1))
        ds_disint = create_dataset(fid, "disintegration", datatype(Float64),
            dataspace(input.box.grid_size..., n_facies, input.time.steps),
            chunk=(input.box.grid_size..., n_facies, 1))
        ds_sedim = create_dataset(fid, "sedim", datatype(Float64),
            dataspace(input.box.grid_size..., n_facies, input.time.steps),
            chunk=(input.box.grid_size..., n_facies, 1))

        results = run_model(input)
        for (step, frame) in enumerate(results)
            ds_prod[:, :, :, step] = frame.production |> in_units_of(u"m")
            ds_disint[:, :, :, step] = frame.disintegration |> in_units_of(u"m")
            ds_sedim[:, :, :, step] = frame.sedimentation |> in_units_of(u"m")
        end
    end
end
```
