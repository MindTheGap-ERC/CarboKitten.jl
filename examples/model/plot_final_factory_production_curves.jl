# plot_final_factory_production_curves.jl
#
# Purpose:
#   Plot production-rate curves for the final 7 CarboKitten factories
#   from the HDF5 production tables saved in data/output/eclepens_withoutca.h5.
#
# Outputs:
#   data/output/production_curves/production_curves_by_factory_stage_midpoints.png
#   data/output/production_curves/production_curves_stage_O.png
#   data/output/production_curves/production_curves_stage_R1.png
#   data/output/production_curves/production_curves_stage_R2.png
#   data/output/production_curves/production_curves_stage_R3.png
#   data/output/production_curves/production_curves_stage_R4.png
#   data/output/production_curves/production_curves_stage_R5.png
#
# Usage from repository root:
#   julia --project=. examples/model/plot_final_factory_production_curves.jl
#
# If CairoMakie is missing:
#   julia --project=.
#   ] add CairoMakie Colors HDF5

using HDF5
using CairoMakie
using Colors
using Printf

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

const H5_FILE = "data/output/eclepens_withoutca.h5"
const OUT_DIR = "data/output/production_curves"

# Final calibrated factory order.
const FACTORY_NAMES = [
    "Ooids",
    "Corals",
    "Mud",
    "Peloids",
    "Bioclasts / intraclasts",
    "Oncoids",
    "Evaporitic mud",
]

const FACTORY_SHORT = [
    "ooids",
    "corals",
    "mud",
    "peloids",
    "bioclasts",
    "oncoids",
    "evap_mud",
]

# Display colours only. They do not affect the model.
const FACTORY_COLORS = [
    colorant"#ffd92f",  # ooids
    colorant"#fb8072",  # corals
    colorant"#6b4f2a",  # mud
    colorant"#b3de69",  # peloids
    colorant"#fdb462",  # bioclasts
    colorant"#80b1d3",  # oncoids
    colorant"#bc80bd",  # evaporitic mud
]

# Final model stage midpoints in elapsed Myr after t0.
# These are used to select representative columns in the production table.
const MODEL_T0_MA = -157.3
const MODEL_DURATION_MYR = 12.3

const STAGE_MIDPOINTS = [
    ("O",  0.5 * (0.00 + 6.40)),
    ("R1", 0.5 * (6.40 + 7.31)),
    ("R2", 0.5 * (7.31 + 9.56)),
    ("R3", 0.5 * (9.56 + 10.75)),
    ("R4", 0.5 * (10.75 + 11.49)),
    ("R5", 0.5 * (11.49 + 12.30)),
]

const STAGE_COLORS = [
    colorant"#1b9e77",
    colorant"#66a61e",
    colorant"#e6ab02",
    colorant"#d95f02",
    colorant"#7570b3",
    colorant"#e7298a",
]

# Set to nothing to plot the full stored depth axis.
const MAX_DEPTH_M = nothing

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

function read_scalar_attribute(g::HDF5.Group, name::AbstractString)
    attrs = HDF5.attributes(g)
    if haskey(attrs, name)
        v = attrs[name]
        try
            return v[]
        catch
            return read(v)
        end
    end
    return nothing
end

function safe_read_vector(g::HDF5.Group, key::AbstractString)
    haskey(g, key) || return nothing
    return vec(Float64.(read(g[key])))
end

function fallback_time_axis(g::HDF5.Group, n_t::Int)
    # If input/t exists but has a different length, use its range.
    # Otherwise use elapsed time from 0 to MODEL_DURATION_MYR.
    for key in ("t", "time", "time_axis", "production_time_axis")
        arr = safe_read_vector(g, key)
        arr === nothing && continue
        if length(arr) == n_t
            return arr
        elseif length(arr) > 1
            return collect(range(minimum(arr), maximum(arr); length=n_t))
        end
    end

    return collect(range(0.0, MODEL_DURATION_MYR; length=n_t))
end

function read_time_axis(input_g::HDF5.Group, facies_g::HDF5.Group, n_t::Int)
    # Prefer time axes stored directly under the facies group.
    for key in ("production_time_axis", "time_axis", "t_axis", "t")
        arr = safe_read_vector(facies_g, key)
        arr === nothing && continue
        if length(arr) == n_t
            return arr
        elseif length(arr) > 1
            return collect(range(minimum(arr), maximum(arr); length=n_t))
        end
    end

    # Then try the input group.
    return fallback_time_axis(input_g, n_t)
end

function stage_target_time(axis::AbstractVector{<:Real}, elapsed_myr::Real)
    # If the stored time axis is absolute geological age, use MODEL_T0 + elapsed.
    # If the stored time axis is elapsed model time, use elapsed directly.
    if minimum(axis) < -1.0 && maximum(axis) < 1.0
        return MODEL_T0_MA + elapsed_myr
    else
        return elapsed_myr
    end
end

function nearest_index(axis::AbstractVector{<:Real}, target::Real)
    _, idx = findmin(abs.(axis .- target))
    return idx
end

function read_production_tables(filename::AbstractString)
    isfile(filename) || error("HDF5 file not found: $filename")

    h5open(filename, "r") do fid
        haskey(fid, "input") || error("HDF5 file has no input group: $filename")
        input_g = fid["input"]

        attr_n = read_scalar_attribute(input_g, "n_facies")
        n_facies = attr_n === nothing ? length(FACTORY_NAMES) : Int(attr_n)

        data = NamedTuple[]

        for i in 1:n_facies
            key = "facies$(i)"
            haskey(input_g, key) || begin
                @warn "Missing input/$key, skipping"
                continue
            end

            fg = input_g[key]
            haskey(fg, "production_table") || begin
                @warn "Missing input/$key/production_table, skipping"
                continue
            end
            haskey(fg, "production_depth_axis") || error("Missing input/$key/production_depth_axis")

            table = Float64.(read(fg["production_table"]))
            depth = vec(Float64.(read(fg["production_depth_axis"])))

            # Expected table order is (n_depth, n_time).
            if size(table, 1) == length(depth)
                table_nd_nt = table
            elseif size(table, 2) == length(depth)
                table_nd_nt = permutedims(table, (2, 1))
            else
                error("For $key, production_table size $(size(table)) is incompatible with depth axis length $(length(depth))")
            end

            n_d, n_t = size(table_nd_nt)
            taxis = read_time_axis(input_g, fg, n_t)

            if MAX_DEPTH_M !== nothing
                keep = findall(d -> d <= MAX_DEPTH_M, depth)
                depth_plot = depth[keep]
                table_plot = table_nd_nt[keep, :]
            else
                depth_plot = depth
                table_plot = table_nd_nt
            end

            name = i <= length(FACTORY_NAMES) ? FACTORY_NAMES[i] : "Factory $(i)"
            short = i <= length(FACTORY_SHORT) ? FACTORY_SHORT[i] : "factory_$(i)"
            color = i <= length(FACTORY_COLORS) ? FACTORY_COLORS[i] : colorant"#444444"

            push!(data, (
                id=i,
                name=name,
                short=short,
                color=color,
                depth=depth_plot,
                table=table_plot,
                time_axis=taxis,
            ))
        end

        return data
    end
end

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

function plot_by_factory_stage_midpoints(data)
    n = length(data)
    ncol = 4
    nrow = ceil(Int, n / ncol)

    fig = Figure(size=(1550, 880), fontsize=16)

    for (idx, d) in enumerate(data)
        row = div(idx - 1, ncol) + 1
        col = mod(idx - 1, ncol) + 1

        ax = Axis(fig[row, col],
            title="$(d.id). $(d.name)",
            xlabel="production [m/Myr]",
            ylabel="depth [m]")

        ax.yreversed = true

        for (k, (stage, elapsed)) in enumerate(STAGE_MIDPOINTS)
            target = stage_target_time(d.time_axis, elapsed)
            ti = nearest_index(d.time_axis, target)
            lines!(ax, d.table[:, ti], d.depth;
                color=STAGE_COLORS[k],
                linewidth=3,
                label=stage)
        end
    end

    # One shared legend.
    elems = [LineElement(color=STAGE_COLORS[k], linewidth=4) for k in eachindex(STAGE_MIDPOINTS)]
    labels = [stage for (stage, elapsed) in STAGE_MIDPOINTS]
    Legend(fig[nrow + 1, 1:ncol], elems, labels, "Stage midpoint",
        orientation=:horizontal, framevisible=false)

    mkpath(OUT_DIR)
    outfile = joinpath(OUT_DIR, "production_curves_by_factory_stage_midpoints.png")
    save(outfile, fig)
    println("Wrote $outfile")
end

function plot_all_factories_for_stage(data, stage_name::AbstractString, elapsed::Real)
    fig = Figure(size=(900, 760), fontsize=18)
    ax = Axis(fig[1, 1],
        xlabel="production [m/Myr]",
        ylabel="depth [m]",
        title="Production curves - stage $(stage_name)")

    ax.yreversed = true

    for d in data
        target = stage_target_time(d.time_axis, elapsed)
        ti = nearest_index(d.time_axis, target)

        lines!(ax, d.table[:, ti], d.depth;
            color=d.color,
            linewidth=3,
            label="$(d.id). $(d.name)")
    end

    Legend(fig[1, 2], ax, "Factory", framevisible=false)

    mkpath(OUT_DIR)
    outfile = joinpath(OUT_DIR, "production_curves_stage_$(stage_name).png")
    save(outfile, fig)
    println("Wrote $outfile")
end

function print_selected_indices(data)
    println("Representative production-table columns:")
    for d in data
        println("  $(d.id). $(d.name)")
        for (stage, elapsed) in STAGE_MIDPOINTS
            target = stage_target_time(d.time_axis, elapsed)
            ti = nearest_index(d.time_axis, target)
            @printf("    %-2s -> column %-4d  stored_time = %.4f\n", stage, ti, d.time_axis[ti])
        end
    end
end

function main()
    data = read_production_tables(H5_FILE)

    isempty(data) && error("No production tables were found in $H5_FILE")

    println("Read $(length(data)) production factories from $H5_FILE")
    for d in data
        println("  $(d.id). $(d.name): table=$(size(d.table)), depth=$(length(d.depth)), time=$(length(d.time_axis))")
    end

    print_selected_indices(data)

    plot_by_factory_stage_midpoints(data)

    for (stage, elapsed) in STAGE_MIDPOINTS
        plot_all_factories_for_stage(data, stage, elapsed)
    end

    println("Done.")
end

main()
