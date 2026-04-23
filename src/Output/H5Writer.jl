module H5Writer

using HDF5
using Unitful
using Unitful: ustrip
using StatsBase: mean

import ...CarboKitten
import ...CarboKitten: Model, AbstractOutput, AbstractInput, OutputSpec, AbstractState, run_model
using ...CarboKitten: time_axis, box_axes
using ...Components.WaterDepth: initial_topography, subsidence_rate_matrix
using ...Components.Production
using ...Utility: in_units_of

using ..Abstract
import ..Abstract: add_data_set, set_attribute, frame_writer, state_writer
export EnvironmentRule,
       build_environment_grid,
       encode_environments,
       env_selector

# --------------------------------------------------
# ENVIRONMENT RULES
# --------------------------------------------------

struct EnvironmentRule
    name::String
    ranges::Dict{Int, Tuple{Float64, Float64}}
    wdepth_range::Union{Nothing, Tuple{Float64,Float64}}
    energy_range::Union{Nothing, Tuple{Float64,Float64}}
end

const env_rules = [
    EnvironmentRule(
        "Lagoon",
        Dict(
            3 => (40.0, 100.0),
            5 => (0.0, 30.0)
        ),
        (0.0, 20.0),
        nothing
    ),

    EnvironmentRule(
        "Shoal",
        Dict(
            5 => (50.0, 100.0),
            3 => (0.0, 40.0)
        ),
        (0.0, 10.0),
        (2.0, 1.0e9)
    ),

    EnvironmentRule(
        "Reef",
        Dict(
            2 => (40.0, 100.0)
        ),
        (0.0, 30.0),
        nothing
    ),

    EnvironmentRule(
        "Deep",
        Dict(),
        (50.0, 5000.0),
        nothing
    )
]

# --------------------------------------------------
# FACIES %
# --------------------------------------------------

function facies_percent_from_ids(fac_ids::AbstractArray{<:Integer,3}; n_facies::Int)
    tot = count(!=(0), fac_ids)
    tot == 0 && return Dict(fid => 0.0 for fid in 1:n_facies)

    return Dict(
        fid => 100.0 * count(==(fid), fac_ids) / tot
        for fid in 1:n_facies
    )
end

function facies_thickness_from_layers(
    facies_layers::AbstractVector,
    thickness_layers::AbstractVector;
    n_facies::Int
)
    props = zeros(Float64, n_facies)
    total = 0.0

    @inbounds for k in eachindex(facies_layers)
        fac = facies_layers[k]
        thick = thickness_layers[k]

        for i in eachindex(fac)
            f = fac[i]
            if f != 0
                h = Float64(thick[i])
                props[f] += h
                total += h
            end
        end
    end

    return Dict(fid => props[fid] for fid in 1:n_facies), total
end

function facies_percent_from_layers(
    facies_layers::AbstractVector,
    thickness_layers::AbstractVector;
    n_facies::Int
)
    props, total = facies_thickness_from_layers(
        facies_layers,
        thickness_layers;
        n_facies = n_facies,
    )

    total == 0.0 && return Dict(fid => 0.0 for fid in 1:n_facies)

    return Dict(fid => 100.0 * props[fid] / total for fid in 1:n_facies)
end

# --------------------------------------------------
# CLASSIFICATION
# --------------------------------------------------

function classify_block(
    props::Dict{Int,Float64},
    rules::Vector{EnvironmentRule};
    wdepth::Float64 = NaN,
    energy::Float64 = NaN
)
    for rule in rules
        ok = true
        for (fid, (lo, hi)) in rule.ranges
            p = get(props, fid, 0.0)
            if p < lo || p > hi
                ok = false
                break
            end
        end
        ok || continue

        if rule.wdepth_range !== nothing
            isnan(wdepth) && continue
            lo, hi = rule.wdepth_range
            (wdepth < lo || wdepth > hi) && continue
        end

        if rule.energy_range !== nothing
            isnan(energy) && continue
            lo, hi = rule.energy_range
            (energy < lo || energy > hi) && continue
        end

        return rule.name
    end

    return "Unclassified"
end

function env_selector(strength, env_rules, env_names;
                      wdepth::Float64 = NaN,
                      energy::Float64 = NaN,
                      nondep_eps = 1e-12)

    total = sum(ustrip.(strength))
    total ≤ nondep_eps && return 0

    props = Dict(
        fid => 100.0 * ustrip(strength[fid]) / total
        for fid in eachindex(strength)
    )

    env_name = classify_block(props, env_rules; wdepth=wdepth, energy=energy)

    if env_name == "Unclassified"
        return length(env_names) - 1
    end

    id = findfirst(==(env_name), env_names)
    return id === nothing ? length(env_names) - 1 : id - 1
end

# --------------------------------------------------
# ENV GRID (TEMP VOXEL -> ENVIRONMENT)
# --------------------------------------------------

function build_environment_grid(
    facies_layers,
    thickness_layers,
    wdepth_layers,
    energy_layers;
    Δi,
    Δj,
    Δk,
    env_rules,
    n_facies
)
    nl = length(facies_layers)
    nx, ny = size(facies_layers[1])

    bi = ceil(Int, nx / Δi)
    bj = ceil(Int, ny / Δj)
    bk = ceil(Int, nl / Δk)

    env_grid = Array{String,3}(undef, bi, bj, bk)

    for I in 1:bi, J in 1:bj, K in 1:bk
        i1 = (I - 1) * Δi + 1
        i2 = min(I * Δi, nx)
        j1 = (J - 1) * Δj + 1
        j2 = min(J * Δj, ny)
        k1 = (K - 1) * Δk + 1
        k2 = min(K * Δk, nl)

        props = Dict(fid => 0.0 for fid in 1:n_facies)
        total_thickness = 0.0
        wd_sum = 0.0
        en_sum = 0.0

        for k in k1:k2
            fac   = facies_layers[k]
            thick = thickness_layers[k]
            wd    = wdepth_layers[k]
            en    = energy_layers[k]

            for ii in i1:i2, jj in j1:j2
                f = fac[ii, jj]
                h = thick[ii, jj]

                if f != 0 && h > 0
                    hh = Float64(h)

                    props[f] += hh
                    total_thickness += hh

                    wd_sum += Float64(wd[ii, jj]) * hh
                    en_sum += Float64(en[ii, jj]) * hh
                end
            end
        end

        if total_thickness == 0
            env_grid[I, J, K] = "Unclassified"
            continue
        end

        for fid in keys(props)
            props[fid] = 100.0 * props[fid] / total_thickness
        end

        mean_wd = wd_sum / total_thickness
        mean_en = en_sum / total_thickness

        env_grid[I, J, K] = classify_block(
            props,
            env_rules;
            wdepth = mean_wd,
            energy = mean_en
        )
    end

    return env_grid
end
#------------------------
# ENCODING
# --------------------------------------------------

function encode_environments(env_grid, env_names)
    map = Dict(name => i for (i, name) in enumerate(env_names))
    encoded = similar(env_grid, Int)

    for i in eachindex(env_grid)
        encoded[i] = get(map, env_grid[i], length(env_names))
    end

    return encoded
end



# --------------------------------------------------
# WRITE ENVIRONMENT
# --------------------------------------------------

function write_environment_block!(fid, encoded, names; Δi, Δj, Δk)
    grp = haskey(fid["input"], "environment_block") ?
          fid["input"]["environment_block"] :
          create_group(fid["input"], "environment_block")

    grp["environment_id"] = encoded

    attrs = attributes(grp)
    attrs["environment_names"] = join(names, ",")
    attrs["Δi"] = Δi
    attrs["Δj"] = Δj
    attrs["Δk"] = Δk
end

# --------------------------------------------------
# WRITE LAYER ARCHIVE (PRIMARY)
# --------------------------------------------------

function write_layer_archive!(fid::HDF5.File, state)
    nl = length(state.layer_facies_hist)
    nl == 0 && return nothing

    parent = fid["input"]
    grp = haskey(parent, "layers") ? parent["layers"] : create_group(parent, "layers")

    nx, ny = size(state.layer_facies_hist[1])

    facies = Array{UInt8,3}(undef, nx, ny, nl)
    thick  = Array{Float32,3}(undef, nx, ny, nl)
    wdepth = Array{Float32,3}(undef, nx, ny, nl)
    energy = Array{Float32,3}(undef, nx, ny, nl)

    for k in 1:nl
        facies[:, :, k] = state.layer_facies_hist[k]
        thick[:, :, k]  = state.layer_thickness_hist[k]
        wdepth[:, :, k] = state.layer_wdepth_hist[k]
        energy[:, :, k] = state.layer_energy_hist[k]
    end

    grp["facies"] = facies
    grp["thickness_m"] = thick
    grp["water_depth_m"] = wdepth
    grp["wave_energy_kwpm"] = energy

    return nothing
end

# --------------------------------------------------
# OUTPUT STRUCT + HEADER
# --------------------------------------------------

mutable struct H5Output <: AbstractOutput
    header
    fid::HDF5.File
end

function make_header(input::AbstractInput)
    t_axis = time_axis(input.time)
    x_axis, y_axis = box_axes(input.box)
    axes = Axes(x=x_axis, y=y_axis, t=t_axis)

    h0 = initial_topography(input)
    sl = input.sea_level.(t_axis)

    return Header(
        tag=input.tag,
        axes=axes,
        Δt=input.time.Δt,
        time_steps=input.time.steps,
        grid_size=input.box.grid_size,
        n_facies=length(input.facies),
        initial_topography=h0,
        sea_level=sl,
        subsidence_rate = subsidence_rate_matrix(input),
        data_sets=Dict(),
        attributes=Dict()
    )
end
function H5Output(input::AbstractInput, filename::String)
    header = make_header(input)
    fid = h5open(filename, "w")
    create_group(fid, "input")

   out = H5Output(header, fid)

    finalizer(out) do x
        close(x.fid)
    end

    return out
end
function H5Output(f, input::AbstractInput, filename::String)
    header = make_header(input)

    h5open(filename, "w") do fid
        create_group(fid, "input")

        out = H5Output(header, fid)
        f(out)
    end
end

# --------------------------------------------------
# HELPERS
# --------------------------------------------------

axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)

slice_str(::Colon) = ":"
slice_str(r::AbstractRange{Int}) = "$(r.start):$(r.stop)"
slice_str(i::Int) = "$(i)"

function get_group(fid::HDF5.File, name::String)
    path = split(name, "/")
    gid = fid["input"]

    for n in path[1:end-1]
        if !haskey(gid, n)
            create_group(gid, n)
        end
        gid = gid[n]
    end

    return gid
end

# --------------------------------------------------
# DATASET ALLOCATION
# --------------------------------------------------

function add_data_set(out::H5Output, name::Symbol, spec::OutputSpec)
    header = out.header
    fid = out.fid

    nf = header.n_facies
    nw = div(header.time_steps, spec.write_interval)
    dims = axis_size.(spec.slice, header.grid_size)

    grp = HDF5.create_group(fid, string(name))
    attrs = attributes(grp)
    attrs["slice"] = join(slice_str.(spec.slice), ",")
    attrs["write_interval"] = spec.write_interval

    HDF5.create_dataset(grp, "production", datatype(Float64),
        dataspace(nf, dims..., nw),
        chunk=(nf, dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "disintegration", datatype(Float64),
        dataspace(nf, dims..., nw),
        chunk=(nf, dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "deposition", datatype(Float64),
        dataspace(nf, dims..., nw),
        chunk=(nf, dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "sediment_thickness", datatype(Float64),
        dataspace(dims..., nw + 1),
        chunk=(dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "compacted_thickness", datatype(Float64),
        dataspace(dims..., nw + 1),
        chunk=(dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "cumulative_subsidence", datatype(Float64),
        dataspace(dims..., nw + 1),
        chunk=(dims..., 1), deflate=3)

    HDF5.create_dataset(grp, "wave_energy", datatype(Float64),
        dataspace(dims..., nw),
        chunk=(dims..., 1), deflate=3)
end

# --------------------------------------------------
# ATTRIBUTE WRITERS
# --------------------------------------------------

function set_attribute(out::H5Output, name::String, value::AbstractArray)
    gid = get_group(out.fid, name)
    tag = split(name, "/")[end]

    if haskey(gid, tag)
        delete_object(gid, tag)
    end

    gid[tag] = value
end

function set_attribute(out::H5Output, name::String, value)
    gid = get_group(out.fid, name)
    attr = attributes(gid)
    tag = split(name, "/")[end]
    attr[tag] = value
end

# --------------------------------------------------
# STATE WRITER
# --------------------------------------------------
reshape_or_fill(x, dims) = x isa AbstractArray ? reshape(x, dims) : fill(x, dims...)

function state_writer(input::AbstractInput, out::H5Output)
    output_spec = input.output
    fid = out.fid
    grid_size = out.header.grid_size

    function (idx::Int, state::AbstractState)
        for (k, v) in output_spec
            dims = axis_size.(v.slice, grid_size)

            if mod(idx - 1, v.write_interval) == 0
                write_id = div(idx - 1, v.write_interval) + 1
                grp = fid[string(k)]

                if write_id <= size(grp["sediment_thickness"], 3)
                    grp["sediment_thickness"][:, :, write_id] =
    reshape_or_fill(state.sediment_height[v.slice...], dims) |> in_units_of(u"m")
                end

                if write_id <= size(grp["compacted_thickness"], 3)
                    grp["compacted_thickness"][:, :, write_id] =
    reshape_or_fill(state.sediment_height[v.slice...], dims) |> in_units_of(u"m")
                end

                if write_id <= size(grp["cumulative_subsidence"], 3)
                    grp["cumulative_subsidence"][:, :, write_id] =
    reshape_or_fill(state.cumulative_subsidence[v.slice...], dims) |> in_units_of(u"m")
                end

                if write_id <= size(grp["wave_energy"], 3) && !isempty(state.energy_hist)
                    hist_idx = min((write_id - 1) * v.write_interval + 1, length(state.energy_hist))
                    grp["wave_energy"][:, :, write_id] =
                        Float64.(state.energy_hist[hist_idx][v.slice...])
                end
            end
        end
    end
end

# --------------------------------------------------
# FRAME WRITER
# --------------------------------------------------

function frame_writer(input::AbstractInput, out::H5Output)
    n_f = out.header.n_facies
    grid_size = out.header.grid_size

    function (idx::Int, frame)
        try_write(tgt, ::Nothing, v) = ()

        function try_write(tgt, src::AbstractArray, v)
            dims = axis_size.(v.slice, grid_size)
            tgt[:, :, :, div(idx - 1, v.write_interval) + 1] +=
                (reshape(src[:, v.slice...], (n_f, dims...)) |> in_units_of(u"m"))
        end

        for (k, v) in input.output
            n_writes = div(input.time.steps, v.write_interval)

            if div(idx - 1, v.write_interval) + 1 <= n_writes
                grp = out.fid[string(k)]
                try_write(grp["production"], frame.production, v)
                try_write(grp["disintegration"], frame.disintegration, v)
                try_write(grp["deposition"], frame.deposition, v)
            end
        end
    end
end

# --------------------------------------------------
# RUN MODEL
# --------------------------------------------------

function run_model(::Type{Model{M}}, input, filename; env_names = String[]) where {M}
    final_state = nothing

    H5Output(input, filename) do output
        final_state = CarboKitten.run_model(Model{M}, input, output)

        # ------------------------------------------
        # Rebuild derived compacted layer products
        # from the true burial / compaction history
        # ------------------------------------------
        if hasproperty(final_state, :compaction) &&
           hasproperty(final_state, :wdepth_hist) &&
           hasproperty(final_state, :energy_hist) &&
           hasproperty(final_state, :block_cube) &&
           hasproperty(final_state, :block_wdepth) &&
           hasproperty(final_state, :block_energy) &&
           hasproperty(final_state, :block_topk)

            Production.extract_compacted_layer_maps!(
                final_state,
                final_state.wdepth_hist,
                final_state.energy_hist
            )

            final_state.block_cube,
            final_state.block_wdepth,
            final_state.block_energy,
            final_state.block_topk =
                Production.project_compacted_layers_to_cube!(
                    final_state,
                    input.depositional_resolution,
                    final_state.wdepth_hist,
                    final_state.energy_hist,
                    final_state.block_cube,
                    final_state.block_wdepth,
                    final_state.block_energy,
                    final_state.block_topk,
                )
        end

        # ------------------------------
        # PRIMARY LAYER ARCHIVE
        # ------------------------------
        if hasproperty(final_state, :layer_facies_hist) &&
           hasproperty(final_state, :layer_thickness_hist) &&
           hasproperty(final_state, :layer_wdepth_hist) &&
           hasproperty(final_state, :layer_energy_hist) &&
           !isempty(final_state.layer_facies_hist)

            write_layer_archive!(output.fid, final_state)

            env_grid = build_environment_grid(
                final_state.layer_facies_hist,
                final_state.layer_thickness_hist,
                final_state.layer_wdepth_hist,
                final_state.layer_energy_hist;
                Δi = 10,
                Δj = 10,
                Δk = 5,
                env_rules = env_rules,
                n_facies = length(input.facies),
            )

            encoded = encode_environments(env_grid, env_names)

            write_environment_block!(
                output.fid,
                encoded,
                env_names;
                Δi = 10,
                Δj = 10,
                Δk = 5,
            )
        end
    end

    return final_state
end
end

