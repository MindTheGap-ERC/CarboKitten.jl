# ~/~ begin <<docs/src/model-alcap.md#src/Model/ALCAPS.jl>>[init]
module ALCAPS

using Unitful
using HDF5
using ProgressBars

using ...BoundaryTrait: Shelf
using ...Config: Box, TimeProperties, axes
using ...SedimentStack: push_sediment!, pop_sediment!
using ...Burgess2013.CA: step_ca
using ...Burgess2013.Production: production_rate
using ...Transport.ActiveLayer: pde_stencil, Rate, Amount
using ...Utility: in_units_of

# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[init]
const Myr = u"Myr"
const m = u"m"

# ~/~ begin <<docs/src/model-alcap.md#alcaps-facies>>[init]
"""
    Facies

Input structure for facies properties.

## Fields

- `viability_range::Tuple{Int,Int}`, range over which cells in the CA stay alive
- `activation_range::Tuple{Int,Int}`, range over which cells in the CA become active if dead
- `maximum_growth_rate`, `extinction_coefficient`, `saturation_intensity`, parameters to the production model [Bosscher1992](@cite)
- `diffusion_coefficient`,  diffusion coefficient in the transport model
"""
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

const FACIES = [
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
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps-sealevel>>[init]
const PERIOD = 0.2Myr
const AMPLITUDE = 4.0m
# ~/~ end

@kwdef struct Input
    tag::String           = "ALCAPS default"
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[init]
    box::Box              = Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m)
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[1]
    time::TimeProperties  = TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1)
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[2]
    ca_interval::Int      = 1
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[3]
    bedrock_elevation     = (x, y) -> -x / 300.0  # (m, m) -> m
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[4]
    sea_level             = t -> AMPLITUDE * sin(2π * t / PERIOD)
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[5]
    subsidence_rate::Rate = 50.0m/Myr
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[6]
    disintegration_rate::Rate = 500.0m/Myr   # same as maximum production rate
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[7]
    insolation::typeof(1.0u"W/m^2") = 400.0u"W/m^2"
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[8]
    sediment_buffer_size::Int     = 50
    # ~/~ end
    # ~/~ begin <<docs/src/model-alcap.md#alcaps-input>>[9]
    depositional_resolution::Amount = 0.5m
    # ~/~ end
    facies::Vector{Facies}    = FACIES
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[1]
"""
    State

## Members

- `time`, absolute time of simulation.
- `ca`, state of the celular automaton.
- `ca_priority`, rotation of activation priority in ca.
- `sediment_height`, the height of the sediment
- `sediment_buffer`, facies composition of sediment

The `sediment_height` ``\\sum_f \\eta`` relates to the water depth ``w`` as follows:

``w = - \\eta_0 + \\sigma t + s(t) - \\sum_f \\eta.``
"""
mutable struct State
    time::typeof(1.0u"Myr")

    ca::Matrix{Int}
    ca_priority::Vector{Int}

    sediment_height::Array{Amount,2}   # x, y
    # sediment_buffer stores fractions, so no units
    sediment_buffer::Array{Float64,4}  # z, facies, x, y
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[2]
"""
    initial_state(input::Input) -> State

Generate the initial state for the model, given the `input`. Returns a `State`.
"""
function initial_state(input)
    sediment_height = zeros(Float64, input.box.grid_size...) * u"m"
    n_facies = length(input.facies)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies, input.box.grid_size...)

    ca = rand(0:n_facies, input.box.grid_size...)
    state = State(0.0u"Myr", ca, 1:n_facies, sediment_height, sediment_buffer)

    for _ = 1:20
        step_ca(input.box, input.facies)(state)
    end

    return state
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[3]
"""
    disintegration(input) -> f!

Prepares the disintegration step. Returns a function `f!(state::State)`. The returned function
modifies the state, popping sediment from the `sediment_buffer` and returns an array of `Amount`.
"""
function disintegration(input)
    n_facies = length(input.facies)
    max_h = input.disintegration_rate * input.time.Δt
    output = Array{Float64, 3}(undef, n_facies, input.box.grid_size...)

    return function(state)
        h = min.(max_h, state.sediment_height)
        state.sediment_height .-= h
        pop_sediment!(state.sediment_buffer, h ./ input.depositional_resolution .|> NoUnits, output)
        return output .* input.depositional_resolution
    end
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[4]
"""
    production(input::Input) -> f

Prepares the production step. Returns a function `f(state::State)`. The returned function
computes the sediment production in an array of `Amount`.
"""
function production(input)
    n_facies = length(input.facies)
    x, y = axes(input.box)
    p(f, w) = production_rate(input.insolation, input.facies[f], w) .* input.time.Δt

    w0 = .- input.bedrock_elevation.(x, y')
    output = Array{Amount, 3}(undef, n_facies, input.box.grid_size...)

    return function(state)
        water_depth = w0 .+ input.sea_level(state.time) .-
            state.sediment_height .+ input.subsidence_rate * state.time
        for f = 1:n_facies
            output[f, :, :] = ifelse.(state.ca .== f, p.(f, water_depth), 0.0m)
        end
        return output
    end
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[5]
"""
    transportation(input::Input) -> f

Prepares the transportation step. Returns a function `f(state::State, active_layer)`,
transporting the active layer, returning a transported `Amount` of sediment.
"""
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
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[6]
"""
    ModelFrame

Output frame of a single iteration of the ALCAPS model.

## Members

- `disintegration`: amount of disintegrated meterial per each facies.
- `production`: amount of produced material per facies.
- `deposition`: transported material, i.e. `transport(disintegration .+ production)`.
- `sediment_height`: resulting sediment height. This is a convenience item, so you don't
  need to recompute the height in post analysis.
"""
struct ModelFrame
    disintegration::Array{Amount,3}    # facies, x, y
    production::Array{Amount,3}
    deposition::Array{Amount,3}
    sediment_height::Array{Amount,2}
end

"""
    run_model(input) -> Channel{ModelFrame}

Runs the ALCAPS model on a given input. Returns a channel of `ModelFrame`.
"""
function run_model(input)
    state = initial_state(input)
    step_ca! = step_ca(input.box, input.facies)
    disintegrate! = disintegration(input)
    produce = production(input)
    transport = transportation(input)

    Channel{ModelFrame}() do ch
        for i in ProgressBar(1:input.time.steps)
            if mod(i, input.ca_interval) == 0
                step_ca!(state)
            end

            p = produce(state)
            d = disintegrate!(state)

            active_layer = p .+ d
            sediment = transport(state, active_layer)

            push_sediment!(state.sediment_buffer, sediment ./ input.depositional_resolution .|> NoUnits)
            state.sediment_height .+= sum(sediment; dims=1)[1,:,:]
            state.time += input.time.Δt

            put!(ch, ModelFrame(d, p, sediment, state.sediment_height))
        end
    end
end

"""
    main(input::Input, output::String)

Runs the ALCAPS model for `input` writing to an HDF5 file given by the `output` path.
"""
function main(input::Input, output::String)
    x, y = axes(input.box)
    t = (0:input.time.steps) .* input.time.Δt

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x) |> in_units_of(u"m")
        gid["y"] = collect(y) |> in_units_of(u"m")
        gid["t"] = t .|> in_units_of(u"Myr")
        gid["bedrock_elevation"] = input.bedrock_elevation.(x, y') |> in_units_of(u"m")
        gid["sea_level"] = input.sea_level.(t) .|> in_units_of(u"m")

        attr = attributes(gid)
        attr["tag"] = input.tag
        attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
        attr["write_interval"] = input.time.write_interval
        attr["time_steps"] = input.time.steps
        attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")
        attr["n_facies"] = length(input.facies)

        n_facies = length(input.facies)
        ds_prod = create_dataset(fid, "production", datatype(Float64),
            dataspace(n_facies, input.box.grid_size..., input.time.steps),
            chunk=(n_facies, input.box.grid_size..., 1))
        ds_disint = create_dataset(fid, "disintegration", datatype(Float64),
            dataspace(n_facies, input.box.grid_size..., input.time.steps),
            chunk=(n_facies, input.box.grid_size..., 1))
        ds_sedim = create_dataset(fid, "deposition", datatype(Float64),
            dataspace(n_facies, input.box.grid_size..., input.time.steps),
            chunk=(n_facies, input.box.grid_size..., 1))
        ds_height = create_dataset(fid, "sediment_height", datatype(Float64),
            dataspace(input.box.grid_size..., input.time.steps),
            chunk=(input.box.grid_size..., 1))

        results = run_model(input)
        for (step, frame) in enumerate(results)
            ds_prod[:, :, :, step] = frame.production |> in_units_of(u"m")
            ds_disint[:, :, :, step] = frame.disintegration |> in_units_of(u"m")
            ds_sedim[:, :, :, step] = frame.deposition |> in_units_of(u"m")
            ds_height[:, :, step] = frame.sediment_height |> in_units_of(u"m")
        end
    end
end
# ~/~ end

end
# ~/~ end
