import CarboKitten.Visualization: facies_colormap

export facies_colormap

function facies_colormap(n_facies; facies_colors=nothing, include_nodeposit=true)

    base = facies_colors === nothing ? FACIES_COLORS : facies_colors

    # Ensure white exists at position 1
    if base[1] != :white
        base = vcat(:white, base)
    end

    required = n_facies + 1

    if length(base) < required
        error("Need $(required) colors (including white), got $(length(base))")
    end

    if include_nodeposit
        cols = base[1:required]                # [white + facies]
        return cgrad(Makie.to_color.(cols), required; categorical=true)
    else
        cols = base[2:required]                # drop white
        return cgrad(Makie.to_color.(cols), n_facies; categorical=true)
    end
end