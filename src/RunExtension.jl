module RunExtension

using CSV
using DataFrames
using Unitful
using Statistics

using ..CarboKitten: run_model, Model
using ..Models: ALCAP
using ..Output.H5Writer: facies_percent_from_ids, facies_percent_from_layers, facies_thickness_from_layers
using ..PatchStats: patch_stats_3D_from_ids
using ..DepthStats: facies_props_depthbins, facies_props_depthbins_thickness

export run_once, run_single_to_csv

# ----------------------------------------------------------
# Utility
# ----------------------------------------------------------

thickness_from_topk(topk::AbstractMatrix{<:Integer}, Δz_m::Real) =
    (mean(topk) * Δz_m,
     minimum(topk) * Δz_m,
     maximum(topk) * Δz_m)

# ----------------------------------------------------------
# Core execution
# ----------------------------------------------------------


function run_once(input;
    run_id::Int = 1,
    n_facies::Int = 5,
    nbins::Int = 35,
    return_state::Bool = true,
    env_names::Vector{String}
)

    state = run_model(
    Model{ALCAP},
    input,
    "run_$(run_id).h5";
    env_names = env_names
)

    nz = maximum(state.block_topk)
    fac_ids = state.block_cube[:, :, 1:nz]

    facies_thickness, total_facies_thickness = facies_thickness_from_layers(
        state.layer_facies_hist,
        state.layer_thickness_hist;
        n_facies = n_facies,
    )
    facies_perc = Dict(
        fid => total_facies_thickness == 0.0 ? 0.0 : 100.0 * facies_thickness[fid] / total_facies_thickness
        for fid in 1:n_facies
    )

    Δz_m = Float64(ustrip(u"m", input.depositional_resolution))

    thickness_mean,
    thickness_min,
    thickness_max =
        thickness_from_topk(state.block_topk, Δz_m)

    patch_stats = patch_stats_3D_from_ids(
        fac_ids;
        n_facies=n_facies,
        r=1
    )

    depthbin_props,
    depthbin_edges =
        facies_props_depthbins(
            fac_ids;
            bin_m=10.0,
            n_facies=n_facies,
            Δz_m=Δz_m
        )


depthbin_props_layers, depthbin_edges_layers =
    facies_props_depthbins_thickness(
        state.layer_facies_hist,
        state.layer_thickness_hist;
        bin_m=10.0,
        n_facies=n_facies,
        Δz_m=Δz_m
    )
    out = (
        run = run_id,
        thickness = thickness_mean,
        thickness_min = thickness_min,
        thickness_max = thickness_max,
        facies_perc = facies_perc,
        facies_thickness = facies_thickness,
        patch_stats = patch_stats,
        depthbin_props = depthbin_props_layers,
        depthbin_edges = depthbin_edges_layers,
        input = input
    )

    return return_state ? merge(out, (state = state,)) : out
end

# ----------------------------------------------------------
# CSV Export
# ----------------------------------------------------------

function run_single_to_csv(input;
    facies_names,
	env_names,
    out_csv,
    fail_csv,
    run_id::Int = 1,
    nbins::Int = 35,
    n_facies::Int = 5
)

    flatten_facies(d::Dict{Int, <:Real}) =
        (; (Symbol("facies_perc_", fid) =>
            Float64(get(d, fid, 0.0))
            for fid in 1:n_facies)...)

    flatten_facies_thickness(d::Dict{Int, <:Real}) =
        (; (Symbol("facies_thickness_", fid) =>
            Float64(get(d, fid, 0.0))
            for fid in 1:n_facies)...)

    function flatten_patches(stats::Dict{Int, NamedTuple})
        cols = Pair{Symbol, Any}[]
        for fid in 1:n_facies
            s = get(stats, fid,
                (n_patches=0,
                 mean_size=0.0,
                 median_size=0.0,
                 max_size=0,
                 total_voxels=0)
            )
            push!(cols, Symbol("patch_n_patches_", fid)    => Int(s.n_patches))
            push!(cols, Symbol("patch_mean_size_", fid)    => Float64(s.mean_size))
            push!(cols, Symbol("patch_median_size_", fid)  => Float64(s.median_size))
            push!(cols, Symbol("patch_max_size_", fid)     => Int(s.max_size))
            push!(cols, Symbol("patch_total_voxels_", fid) => Int(s.total_voxels))
        end
        return (; cols...)
    end

    function flatten_depthbins(props::Dict{Int, Dict{Int, Float64}})
        cols = Pair{Symbol, Any}[]
        for b in 1:nbins
            d = get(props, b, Dict{Int, Float64}())
            for fid in 1:n_facies
                push!(cols,
                    Symbol("depthbin_", b, "_facies_", fid) =>
                        Float64(get(d, fid, 0.0))
                )
            end
        end
        return (; cols...)
    end

    try
        r = run_once(
    input;
    run_id=run_id,
    n_facies=n_facies,
    nbins=nbins,
    return_state=true,
    env_names=env_names
)

        base = (
            run = r.run,
            thickness     = r.thickness,
            thickness_min = r.thickness_min,
            thickness_max = r.thickness_max
        )

        row = merge(
            base,
            flatten_facies(r.facies_perc),
            flatten_facies_thickness(r.facies_thickness),
            flatten_patches(r.patch_stats),
            flatten_depthbins(r.depthbin_props)
        )

        CSV.write(
            out_csv,
            DataFrame([row]);
            append=isfile(out_csv)
        )

        return (out_csv = out_csv,
                state = r.state,
                input = input)

    catch e
        CSV.write(
            fail_csv,
            DataFrame([(run=run_id,
                        error=sprint(showerror, e))]);
            append=isfile(fail_csv)
        )
        rethrow()
    end
end

end # module