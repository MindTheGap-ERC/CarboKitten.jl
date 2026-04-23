module DepthStats

using Statistics: mean
using StatsBase: countmap   # ← only if StatsBase already in Project.toml

export facies_props_depthbins, facies_props_depthbins_thickness

function facies_props_depthbins(
    fac_ids::AbstractArray{<:Integer,3};
    bin_m::Real = 10.0,
    n_facies::Int = 5,
    Δz_m::Real
)

    nz = size(fac_ids, 3)

    # depth of voxel centers
    z_m = ((1:nz) .- 0.5) .* Δz_m
    zmax = maximum(z_m)

    edges = collect(0.0:bin_m:(ceil(zmax/bin_m)*bin_m))
    nb = length(edges) - 1

    out = Dict{Int, Dict{Int, Float64}}()

    for b in 1:nb
        zlo, zhi = edges[b], edges[b+1]

        zidx = findall(z -> (z ≥ zlo) && (z < zhi), z_m)
        isempty(zidx) && continue

        slab = fac_ids[:, :, zidx]

        tot = count(!=(0), slab)
        cm  = countmap(vec(slab))

        props = Dict{Int, Float64}()

        for fid in 1:n_facies
            props[fid] =
                tot == 0 ? 0.0 :
                100.0 * get(cm, fid, 0) / tot
        end

        out[b] = props
    end

    return out, edges
end

function facies_props_depthbins_thickness(
    facies_layers::AbstractVector{<:AbstractMatrix{<:Integer}},
    thickness_layers::AbstractVector{<:AbstractMatrix{<:Real}};
    bin_m::Real=10.0,
    n_facies::Int=5,
    Δz_m::Real
)
    length(facies_layers) == length(thickness_layers) ||
        throw(ArgumentError("`facies_layers` and `thickness_layers` must have the same length"))

    nlayers = length(facies_layers)
    max_depth = nlayers * Δz_m
    nbins = ceil(Int, max_depth / bin_m)
    edges = collect(0.0:bin_m:(nbins * bin_m))

    props = zeros(Float64, nbins, n_facies)

    for b in 1:nbins
        zmin = edges[b]
        zmax = edges[b + 1]

        totals = zeros(Float64, n_facies)
        total_thickness = 0.0

        for k in 1:nlayers
            zmid = (k - 0.5) * Δz_m
            if !(zmin <= zmid < zmax)
                continue
            end

            fac = facies_layers[k]
            thick = thickness_layers[k]
            size(fac) == size(thick) ||
                throw(ArgumentError("Layer $k has mismatched facies/thickness shapes"))

            for j in axes(fac, 2), i in axes(fac, 1)
                f = Int(fac[i, j])
                h = Float64(thick[i, j])

                if 1 <= f <= n_facies && h > 0
                    totals[f] += h
                    total_thickness += h
                end
            end
        end

        if total_thickness > 0
            props[b, :] .= 100.0 .* totals ./ total_thickness
        end
    end

    return props, edges
end
end # module
