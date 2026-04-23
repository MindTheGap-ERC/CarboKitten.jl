# ~/~ begin <<docs/src/model-alcap.md#src/Models/ALCAP.jl>>[init]
@compose module ALCAP
    @mixin Tag, Output, CAProduction, ActiveLayer, InitialSediment, Diagnostics

    using ..Common
    using ..CAProduction: production
    using ..TimeIntegration
    using ..WaterDepth: water_depth, subsidence_rate_matrix
    using ...Output: Frame
    using ModuleMixins: @for_each
    using Unitful: ustrip
    using Unitful: NoUnits
    using ...WaveField
    using GeometryBasics

    export Input, Facies

    

    function initial_state(input::AbstractInput)
        ca_state = CellularAutomaton.initial_state(input)
        for _ in 1:20
            CellularAutomaton.step!(input)(ca_state)
        end

        sediment_height = zeros(Height, input.box.grid_size...)
        sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
        active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)
        cumulative_subsidence = zeros(Height, input.box.grid_size...)
        cumulative_subsidence_hist = Matrix{Height}[]

        nx, ny = input.box.grid_size
        nz_compaction = n_facies(input) * input.time.steps
        compaction = Production.CompactionState(nx, ny, nz_compaction)

        nz0 = 400

        block_cube = zeros(UInt8, nx, ny, nz0)
        block_topk = zeros(Int, nx, ny)

        production_hist = Vector{Array{Float64, 3}}()
        disintegration_hist = Vector{Array{Float64, 3}}()
        deposition_hist = Vector{Array{Float64, 3}}()
        sediment_thickness_hist = Vector{Array{Float64, 2}}()
        time_hist = Float64[]
        block_wdepth = zeros(Float32, nx, ny, nz0)
        block_energy = zeros(Float32, nx, ny, nz0)
        wdepth_hist = Vector{Matrix{Float64}}()
        energy_hist = Vector{Matrix{Float32}}()
        layer_facies_hist = Vector{Matrix{UInt8}}()
        layer_thickness_hist = Vector{Matrix{Float32}}()
        layer_wdepth_hist = Vector{Matrix{Float32}}()
        layer_energy_hist = Vector{Matrix{Float32}}()

        push!(wdepth_hist, zeros(Float64, input.box.grid_size...))
        push!(energy_hist, zeros(Float32, input.box.grid_size...))

        state = State(
            step = 0,
            sediment_height = sediment_height,
            sediment_buffer = sediment_buffer,
            active_layer = active_layer,
            ca = ca_state.ca,
            ca_priority = ca_state.ca_priority,

            cumulative_subsidence = cumulative_subsidence,
            cumulative_subsidence_hist = cumulative_subsidence_hist,

            block_cube = block_cube,
            block_topk = block_topk,

            production_hist = production_hist,
            disintegration_hist = disintegration_hist,
            deposition_hist = deposition_hist,
            sediment_thickness_hist = sediment_thickness_hist,
            block_wdepth = block_wdepth,
            time_hist = time_hist,
            wdepth_hist = wdepth_hist,
            energy_hist = energy_hist,
            block_energy = block_energy,
            layer_facies_hist = layer_facies_hist,
            layer_thickness_hist = layer_thickness_hist,
            layer_wdepth_hist = layer_wdepth_hist,
            layer_energy_hist = layer_energy_hist,
            compaction = compaction,
        )

        InitialSediment.push_initial_sediment!(input, state)

        return state
    end

    function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    transport! = ActiveLayer.transporter(input)
    local_water_depth = water_depth(input)
    pf = cementation_factor(input)
    dtf = input.disintegration_transfer
    debug = input.diagnostics
    Kcomp = 2

    function (state::State)
        if debug
            @debug "step: " state.step
            @debug "   current active layer ambitus: " extrema(state.active_layer)
        end

        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        x_coords, y_coords = box_axes(input.box)
        current_time = input.time.t0 + state.step * input.time.Δt

        sr = subsidence_rate_matrix(input)
        mult = ones(Float64, size(sr))

        for mod_ in input.subsidence_modifiers
            mult .*= mod_.multiplier.(x_coords, y_coords', current_time)
        end

        increment = sr .* mult .* input.time.Δt
        state.cumulative_subsidence .+= increment
        push!(state.cumulative_subsidence_hist, copy(state.cumulative_subsidence))

        # compaction only to compute runtime water depth
        raw_height = copy(state.sediment_height)

        if state.step > 0 && mod(state.step, Kcomp) == 0
            Production.compact_column!(input, state)
            Production.update_sediment_height_from_compaction!(state)
        end

        wd = local_water_depth(state)

        # restore morphodynamic surface immediately
        state.sediment_height .= raw_height

        push!(state.wdepth_hist, Float64.(ustrip.(wd ./ u"m")))

        if input.wave_model !== nothing
            nx, ny = size(wd)
            energy_kwpm = zeros(Float32, nx, ny)
            x0 = first(x_coords)

            for i in 1:nx, j in 1:ny
                h = max(wd[i, j], 0.0u"m")
                w = x_coords[i] + y_coords[j]
                x_fetch = abs(x_coords[i] - x0)

                (_, _, p) = WaveField.wave_physics_at_cell(
                    input.wave_model, w, current_time, h, x_fetch
                )
                energy_kwpm[i, j] = Float32(ustrip(u"kW/m", p))
            end

            push!(state.energy_hist, energy_kwpm)
        else
            push!(state.energy_hist, zeros(Float32, size(wd)...))
        end

        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)

        if debug
            @debug "   post-production ambitus: " extrema(state.active_layer)
        end

        transport!(state)

        if debug
            @debug "   post-transport ambitus: " extrema(state.active_layer)
        end

        deposit = pf .* state.active_layer
        deposited = permutedims(deposit, (2, 3, 1)) .|> ustrip

        Production.bury_deposition!(state, deposited, state.step)

        push_sediment!(
            state.sediment_buffer,
            deposit ./ input.depositional_resolution .|> NoUnits
        )
state.sediment_height .+= dropdims(sum(deposit; dims=1), dims=1)
        state.active_layer .-= deposit

        push!(state.production_hist, ustrip.(p ./ u"m"))
        push!(state.disintegration_hist, ustrip.(d ./ u"m"))
        push!(state.deposition_hist, ustrip.(deposit ./ u"m"))
        push!(state.sediment_thickness_hist, ustrip.(state.sediment_height ./ u"m"))
        push!(state.time_hist, ustrip(state.step * input.time.Δt / u"Myr"))

        state.step += 1

        return Frame(
            production = p,
            disintegration = d,
            deposition = deposit,
        )
    end
end

    function write_header(input::AbstractInput, output::AbstractOutput)
        @for_each(P -> P.write_header(input, output), PARENTS)
    end
end
# ~/~ end