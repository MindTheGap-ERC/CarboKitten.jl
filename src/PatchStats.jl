module PatchStats

using Statistics
using StatsBase   # if you use countmap
using ImageMorphology  # or whatever provides label_components

export patch_stats_3D_from_ids

function patch_stats_3D_from_ids(fac_ids::AbstractArray{<:Integer,3};
                                 n_facies::Int=5, r::Int=1)

    nx, ny, nz = size(fac_ids)
    stats = Dict{Int, NamedTuple}()

    for fid in 1:n_facies
        mask = (fac_ids .== fid)

        if !any(mask)
            stats[fid] = (n_patches=0,
                          mean_size=0.0,
                          median_size=0.0,
                          max_size=0,
                          total_voxels=0)
            continue
        end

        lab = label_components(mask; dims=3, r=r)
        counts = countmap(lab)
        ps = Int[v for (k,v) in counts if k != 0]

        stats[fid] = (
            n_patches    = length(ps),
            mean_size    = mean(ps),
            median_size  = median(ps),
            max_size     = maximum(ps),
            total_voxels = sum(ps)
        )
    end

    return stats
end

end
