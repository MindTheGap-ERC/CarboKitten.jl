# ~/~ begin <<docs/src/components/production.md#src/Components/Production.jl>>[init]
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using HDF5
using Logging
using Unitful: ustrip


export production_rate, capped_production, uniform_production

@kwdef struct Facies <: AbstractFacies
    maximum_growth_rate::Rate = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::Intensity = 1.0u"W/m^2"
    depth_knots::Vector{Tuple{typeof(1.0u"m"), Float64}} =
        Tuple{typeof(1.0u"m"), Float64}[]

    time_windows::Vector{Tuple{typeof(1.0u"Myr"), typeof(1.0u"Myr"), Float64}} =
        Tuple{typeof(1.0u"Myr"), typeof(1.0u"Myr"), Float64}[]
end


@kwdef struct Input <: AbstractInput
    insolation
end


@kwdef mutable struct CompactionState
    layer_facies::Array{Int16,3}
    layer_thickness0::Array{Float64,3}
    layer_thickness::Array{Float64,3}
    layer_step::Array{Int32,3}
    n_layers::Array{Int16,2}
end

function CompactionState(nx::Int, ny::Int, nzmax::Int)
    CompactionState(
        fill(Int16(0), nx, ny, nzmax),
        fill(0.0, nx, ny, nzmax),
        fill(0.0, nx, ny, nzmax),
        fill(Int32(0), nx, ny, nzmax),   # ← layer_step
        fill(Int16(0), nx, ny),
    )
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
    set_attribute(output, "facies$(i)/maximum_growth_rate",
        f.maximum_growth_rate |> in_units_of(u"m/Myr"))
set_attribute(output, "facies$(i)/extinction_coefficient",
    f.extinction_coefficient |> in_units_of(u"m^-1"))
set_attribute(output, "facies$(i)/saturation_intensity",
    f.saturation_intensity |> in_units_of(u"W/m^2"))
    # depth knots
    if !isempty(f.depth_knots)
        depths = first.(f.depth_knots) |> in_units_of(u"m")
        mults  = last.(f.depth_knots)
        set_attribute(output, "facies$(i)/depth_knots/depth_m", depths)
        set_attribute(output, "facies$(i)/depth_knots/mult", mults)
    end

    # time windows
    if !isempty(f.time_windows)
        t1 = first.(f.time_windows) |> in_units_of(u"Myr")
        t2 = getindex.(f.time_windows, 2) |> in_units_of(u"Myr")
        mm = last.(f.time_windows)
        set_attribute(output, "facies$(i)/time_windows/t1_Myr", t1)
        set_attribute(output, "facies$(i)/time_windows/t2_Myr", t2)
        set_attribute(output, "facies$(i)/time_windows/mult", mm)
    end
end


end

# ~/~ begin <<docs/src/components/production.md#component-production-rate>>[init]

function depth_multiplier(f, w)
    k = f.depth_knots
    isempty(k) && return 1.0
    isempty(k) && return 1.0
    w < k[1][1] && return 0.0
    w > k[end][1] && return 0.0
    for i in 1:length(k)-1
        (w1,m1) = k[i]; (w2,m2) = k[i+1]
        if w1 <= w <= w2
            α = ustrip((w - w1)/(w2 - w1))
            return (1-α)*m1 + α*m2
        end
    end
    return 1.0
end

function time_multiplier(f, t)
    m = 1.0
    for (t1,t2,mm) in f.time_windows
        if t1 <= t <= t2
            m = mm   # last match wins
        end
    end
    return m
end

function porosity_at_depth(f, z)
    ϕ = f.compaction_curve(z)

    if !isfinite(ϕ)
        error("Compaction curve returned non-finite porosity at z = $z")
    end

    if !(0.0 <= ϕ < 1.0)
        error("Compaction curve returned invalid porosity $ϕ at z = $z; expected 0 <= ϕ < 1")
    end

    return ϕ
end


function compacted_thickness(H0, ϕ0, ϕ)
    return H0 * (1.0 - ϕ0) / max(1.0 - ϕ, 1e-8)
end



function compact_column!(input, state)
    comp = state.compaction
    nx, ny = size(comp.n_layers)

    for i in 1:nx, j in 1:ny
        z_top = 0.0

        for k in reverse(1:comp.n_layers[i, j])
            fidx = comp.layer_facies[i, j, k]
            H0   = comp.layer_thickness0[i, j, k]

            facies = input.facies[fidx]
            ϕ0 = facies.depositional_porosity

            z_mid = z_top + 0.5 * H0
            ϕ = porosity_at_depth(facies, z_mid)

            Hc = compacted_thickness(H0, ϕ0, ϕ)
            comp.layer_thickness[i, j, k] = Hc

            z_top += Hc
        end
    end
end


function update_sediment_height_from_compaction!(state)
    comp = state.compaction
    nx, ny = size(comp.n_layers)

    for i in 1:nx, j in 1:ny
        total = 0.0

        for k in 1:comp.n_layers[i, j]
            total += comp.layer_thickness[i, j, k]
        end

        active = sum(state.active_layer[:, i, j])
        state.sediment_height[i, j] = (total + ustrip(active / u"m")) * u"m"
    end
end


function merge_layers!(state)
    comp = state.compaction
    nx, ny = size(comp.n_layers)

    for i in 1:nx, j in 1:ny
        n = comp.n_layers[i, j]
        n <= 1 && continue

        new_k = 1

        for k in 2:n
            f_prev = comp.layer_facies[i, j, new_k]
            f_curr = comp.layer_facies[i, j, k]

            if f_curr == f_prev
                # merge into previous layer
                comp.layer_thickness0[i, j, new_k] += comp.layer_thickness0[i, j, k]
                comp.layer_thickness[i, j, new_k]  += comp.layer_thickness[i, j, k]
            else
                # move layer forward
                new_k += 1
                comp.layer_facies[i, j, new_k] = f_curr
                comp.layer_thickness0[i, j, new_k] = comp.layer_thickness0[i, j, k]
                comp.layer_thickness[i, j, new_k]  = comp.layer_thickness[i, j, k]
            end
        end

        comp.n_layers[i, j] = new_k
    end
end


function extract_compacted_layer_maps!(state, wdepth_hist, energy_hist)
    comp = state.compaction
    nx, ny = size(comp.n_layers)
    max_layers = maximum(comp.n_layers)

    empty!(state.layer_facies_hist)
    empty!(state.layer_thickness_hist)
    empty!(state.layer_wdepth_hist)
    empty!(state.layer_energy_hist)

    for k in 1:max_layers
        facies_map = zeros(UInt8, nx, ny)
        thickness_map = zeros(Float32, nx, ny)
        wdepth_map = zeros(Float32, nx, ny)
        energy_map = zeros(Float32, nx, ny)

        for i in 1:nx, j in 1:ny
            if k <= comp.n_layers[i, j]
                facies_map[i, j] = UInt8(comp.layer_facies[i, j, k])
                thickness_map[i, j] = Float32(comp.layer_thickness[i, j, k])
                s = comp.layer_step[i, j, k]

		hist_idx = clamp(Int(s) + 1, 1,
    		min(length(wdepth_hist), length(energy_hist))
		)

		wdepth_map[i, j] = Float32(wdepth_hist[hist_idx][i, j])
		energy_map[i, j] = Float32(energy_hist[hist_idx][i, j])
            end
        end

        push!(state.layer_facies_hist, facies_map)
        push!(state.layer_thickness_hist, thickness_map)
        push!(state.layer_wdepth_hist, wdepth_map)
        push!(state.layer_energy_hist, energy_map)
    end
end






function project_compacted_layers_to_cube!(
    state,
    Δz,
    wdepth_hist,
    energy_hist,
    block_cube,
    block_wdepth,
    block_energy,
    block_topk,
)
    comp = state.compaction
    nx, ny = size(comp.n_layers)
    nz = size(block_cube, 3)
    Δz_m = ustrip(Δz / u"m")

    # --------------------------------------------------
    # First pass: compute required nz
    # --------------------------------------------------
    needed_nz = 0
    for i in 1:nx, j in 1:ny
        carry = 0.0
        kk = 0

        for k in 1:comp.n_layers[i, j]
            H = comp.layer_thickness[i, j, k]
            h_total = H + carry
            nvox = Int(floor(h_total / Δz_m))
            carry = h_total - nvox * Δz_m
            kk += max(nvox, 0)
        end

        needed_nz = max(needed_nz, kk)
    end

    # --------------------------------------------------
    # Resize if needed
    # --------------------------------------------------
    if needed_nz > nz
        new_nz = max(needed_nz, Int(ceil(1.2 * nz)) + 1)

        block_cube   = zeros(UInt8, nx, ny, new_nz)
        block_wdepth = zeros(Float32, nx, ny, new_nz)
        block_energy = zeros(Float32, nx, ny, new_nz)
    end

    # --------------------------------------------------
    # Reset arrays
    # --------------------------------------------------
    fill!(block_cube, 0x00)
    fill!(block_wdepth, 0.0f0)
    fill!(block_energy, 0.0f0)
    fill!(block_topk, 0)

    # --------------------------------------------------
    # Projection
    # --------------------------------------------------
    for i in 1:nx, j in 1:ny
        carry = 0.0
        kk = 0

        for k in 1:comp.n_layers[i, j]
            f = comp.layer_facies[i, j, k]
            H = comp.layer_thickness[i, j, k]
            s = comp.layer_step[i, j, k]

            h_total = H + carry
            nvox = Int(floor(h_total / Δz_m))
            carry = h_total - nvox * Δz_m

            nvox <= 0 && continue

            hist_idx = clamp(Int(s) + 1, 1, min(length(wdepth_hist), length(energy_hist)))

            wd = Float32(wdepth_hist[hist_idx][i, j])
            en = Float32(energy_hist[hist_idx][i, j])

            for _ in 1:nvox
                kk += 1
                if kk > size(block_cube, 3)
                    error("Projected compacted cube still exceeds allocated nz at ($(i), $(j))")
                end

                block_cube[i, j, kk] = UInt8(f)
                block_wdepth[i, j, kk] = wd
                block_energy[i, j, kk] = en
            end
        end

        block_topk[i, j] = kk
    end

    return block_cube, block_wdepth, block_energy, block_topk
end



function bury_deposition!(state, deposited, step)
    comp = state.compaction
    nx, ny, nf = size(deposited)

    for i in 1:nx, j in 1:ny, f in 1:nf
        H0 = deposited[i, j, f]

        H0 <= 1e-6 && continue

        k = comp.n_layers[i, j] + 1

        if k > size(comp.layer_thickness0, 3)
            error("Maximum number of buried layers exceeded at ($i, $j)")
        end

        comp.layer_facies[i, j, k] = f
        comp.layer_thickness0[i, j, k] = H0
        comp.layer_thickness[i, j, k] = H0
        comp.n_layers[i, j] = k
        comp.layer_step[i, j, k] = step
    end
end




function production_rate(f, insolation, w::Height, t::Time)
    if isempty(f.depth_knots) && isempty(f.time_windows)
        I = insolation / f.saturation_intensity
        x = w * f.extinction_coefficient
        return x > 0.0 ? f.maximum_growth_rate * tanh(I * exp(-x)) : 0.0u"m/Myr"
    end

    return f.maximum_growth_rate * depth_multiplier(f, w) * time_multiplier(f, t)
end

function capped_production(f, insolation, w::Height, t::Time, dt::Time)
    p = production_rate(f, insolation, w, t) * dt
    return min(p, max(0.0u"m", w))
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies
    dt = input.time.Δt
    tprop = input.time

    p(state::AbstractState, wd::AbstractMatrix) = begin
        t = time(tprop, state)
        I = insolation_func(state)
        capped_production.(facies[:, na, na], I, wd[na, :, :], t, dt)
    end

    p(state::AbstractState) = p(state, w(state))
    return p
end


end
# ~/~ end
