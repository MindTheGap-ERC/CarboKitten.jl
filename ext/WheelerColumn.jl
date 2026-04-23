using Makie
using Unitful: ustrip, @u_str
using Statistics
import CarboKitten.Visualization: wheeler_column, facies_colormap
export wheeler_column


function wheeler_column(
    state;
    # spatial selection
    x::Int,
    y::Int = nothing,
    Δi::Int = 1,
    Δj::Int = 1,

    # category definition
    category_names::Vector{String},
    category_colors,
    selector::Function,   # (strength, wdepth, energy_kwpm) -> category id
    nondep_eps = 1e-12,

    title_str::String = "",
    colorbar_label::String = ""
)

    dep_hist = state.deposition_hist
    nt = length(dep_hist)
    nt > 0 || error("state.deposition_hist is empty.")

    nf, nx, ny = size(dep_hist[1])
    1 ≤ x ≤ nx || error("Invalid x index.")

    # ---- define spatial window ----
    if y === nothing
        # column across all y (facies x-profile case)
        i1, i2 = x, x
        j1, j2 = 1, ny
    else
        1 ≤ y ≤ ny || error("Invalid y index.")

        I = ceil(Int, x / Δi)
        J = ceil(Int, y / Δj)

        i1 = (I-1)*Δi + 1
        i2 = min(I*Δi, nx)

        j1 = (J-1)*Δj + 1
        j2 = min(J*Δj, ny)
    end

    # ---- time axis ----
    tmyr = try
        Float64.(ustrip.(u"Myr", state.time_hist))
    catch
        Float64.(state.time_hist)
    end

    # ---- dominant category per timestep ----
    dom_code = zeros(Float32, nt)

    for it in 1:nt
        dep = dep_hist[it]
        block = dep[:, i1:i2, j1:j2]
		wd2d = state.wdepth_hist[it]
		wd_block = wd2d[i1:i2, j1:j2]
		mean_wd = mean(wd_block)
		
		e2d = state.energy_hist[it]
		e_block = e2d[i1:i2, j1:j2]
		mean_e = mean(e_block)

        strength = vec(sum(block; dims=(2,3)))

        if sum(strength) ≤ nondep_eps
            dom_code[it] = 0.0f0
        else
            dom_code[it] = try
				Float32(selector(strength, mean_wd, mean_e))
			catch err
				err isa MethodError || rethrow()
				Float32(selector(strength, mean_wd))  # backward compatible
			end
        end
    end

    # ---- colormap ----
    ncat = length(category_names)

    cmap = facies_colormap(
    ncat;
    facies_colors = category_colors,
    include_nodeposit = true
)

    img = repeat(reshape(dom_code, 1, :), 2, 1)

    fig = Figure(size=(360, 850))
    ax = Axis(fig[1, 1];
        ylabel = "Time (Myr)",
        title  = isempty(title_str) ?
            "Wheeler column" :
            title_str,
        xticksvisible = false,
        xticklabelsvisible = false
    )

    heatmap!(ax, [0.0, 1.0], tmyr, img;
        colormap   = cmap,
        colorrange = (-0.5f0, ncat - 0.5f0)
    )

    Colorbar(fig[1, 2];
    colormap = cmap,
    limits   = (-0.5, ncat + 0.5),
    ticks    = (0:ncat, ["No deposit"; category_names]),
    label    = colorbar_label
)

    return fig
end
