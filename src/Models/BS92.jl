# ~/~ begin <<docs/src/bosscher-1992.md#src/Models/BS92.jl>>[init]
@compose module BS92
@mixin Tag, Diagnostics, Output, Production

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using ...Output: Frame
using Unitful: ustrip, unit

export Input, Facies

function initial_state(input::Input)
    nx, ny = input.box.grid_size

    sediment_height = zeros(Height, nx, ny)
    cumulative_subsidence = zeros(Height, nx, ny)

    cumulative_subsidence_hist = [zeros(Height, nx, ny)]
    wdepth_hist = [zeros(Float64, nx, ny)]
    energy_hist = [zeros(Float32, nx, ny)]

    block_wdepth = zeros(Float32, nx, ny, 1)
    block_energy = zeros(Float32, nx, ny, 1)

    return State(
        step = 0,
        sediment_height = sediment_height,
        cumulative_subsidence = cumulative_subsidence,
        cumulative_subsidence_hist = cumulative_subsidence_hist,
        wdepth_hist = wdepth_hist,
        block_wdepth = block_wdepth,
        energy_hist = energy_hist,
        block_energy = block_energy,
    )
end

function step!(input::Input)
    τ = uniform_production(input)
    function (state::State)
        prod = τ(state)
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
        state.step += 1
        return Frame(
            production=prod,
            deposition=prod)
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
set_attribute(output, "subsidence_rate", ustrip(input.subsidence_rate))
set_attribute(output, "subsidence_rate_units", string(unit(input.subsidence_rate)))
end

end
# ~/~ end
