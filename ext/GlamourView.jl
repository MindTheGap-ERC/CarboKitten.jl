module GlamourView

import CarboKitten.Visualization: glamour_view, glamour_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, DataVolume, datavolume_from_state
import CarboKitten.Components.Production

using Makie
using Unitful

export glamour_view, glamour_view!

# ----------------------------------------------------------------------
# 1. Original method (DataVolume)
# ----------------------------------------------------------------------
function glamour_view!(ax::Makie.Axis3, header::Header, data::DataVolume; colormap=Reverse(:speed))

    x = header.axes.x[data.slice[1]] |> in_units_of(u"km")
    y = header.axes.y[data.slice[2]] |> in_units_of(u"km")
    xy_aspect = x[end] / y[end]

    ax.aspect = (xy_aspect, 1, 0.001)
    ax.azimuth = -π/3

    # --------------------------------------------------
    # PICK DATASET (NEW LOGIC)
    # --------------------------------------------------
    thickness = if hasproperty(data, :compacted_thickness)
        data.compacted_thickness
    else
        data.sediment_thickness
    end

    n_steps = size(thickness, 3)
    grid_size = (length(x), length(y))

    steps_between = 2
    selected_steps = [1, ((1:steps_between) .* n_steps .÷ (steps_between + 1))..., n_steps]

    bedrock = header.initial_topography .-
              (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate

    result = Array{Float64, 3}(undef, grid_size..., length(selected_steps))

    for (i, j) in enumerate(selected_steps)
        result[:, :, i] = (thickness[:,:,j] .+ bedrock) |> in_units_of(u"m")
    end

    # bottom
    surface!(ax, x, y, result[:,:,1];
        color = ones(grid_size),
        colormap = :grays
    )

    # intermediate
    for s in eachslice(result[:,:,2:end-1], dims=3)
        surface!(ax, x, y, s;
            colormap = (colormap, 0.7)
        )
    end

    # top
    surface!(ax, x, y, result[:,:,end];
        colormap = colormap
    )

    # outline
    lines!(ax, x, zeros(grid_size[1]), result[:, 1, end];
        color = (:white, 0.5), linewidth = 1)

    lines!(ax, fill(x[end], grid_size[2]), y, result[end, :, end];
        color = (:white, 0.5), linewidth = 1)
end

# ----------------------------------------------------------------------
# 2. NEW method (compaction-aware, state-based)
# ----------------------------------------------------------------------
function glamour_view!(ax::Makie.Axis3, header::Header, input, state; colormap=Reverse(:speed))
    x = header.axes.x |> in_units_of(u"km")
    y = header.axes.y |> in_units_of(u"km")

    xy_aspect = x[end] / y[end]

    ax.aspect = (xy_aspect, 1, 0.001)
    ax.azimuth = -π/3

    n_steps = state.step
    grid_size = size(state.sediment_height)

    steps_between = 2
    selected_steps = [1, ((1:steps_between) .* n_steps .÷ (steps_between + 1))..., n_steps]

    surfaces = build_selected_surfaces(state, input, selected_steps)

    # bottom surface
    surface!(ax, x, y, surfaces[1];
        color = ones(grid_size),
        colormap = :grays
    )

    # intermediate surfaces
    for s in surfaces[2:end-1]
        surface!(ax, x, y, s;
            colormap = (colormap, 0.7)
        )
    end

    # top surface
    surface!(ax, x, y, surfaces[end];
        colormap = colormap
    )

    # outline
    lines!(ax, x, zeros(grid_size[1]), surfaces[end][:, 1];
        color = (:white, 0.5), linewidth = 1)

    lines!(ax, fill(x[end], grid_size[2]), y, surfaces[end][end, :];
        color = (:white, 0.5), linewidth = 1)
end

# ----------------------------------------------------------------------
# 3. Wrapper (DataVolume)
# ----------------------------------------------------------------------
function glamour_view(header::Header, data::DataVolume; colormap=Reverse(:speed))
    fig = Figure()
    ax = Axis3(fig[1, 1])
    glamour_view!(ax, header, data; colormap=colormap)
    return fig, ax
end

# ----------------------------------------------------------------------
# 4. Wrapper (state + compaction)
# ----------------------------------------------------------------------
function glamour_view(header::Header, input, state; colormap=Reverse(:speed))
    fig = Figure()
    ax = Axis3(fig[1, 1])
    glamour_view!(ax, header, input, state; colormap=colormap)
    return fig, ax
end

# ----------------------------------------------------------------------
# 5. Compaction-aware surface reconstruction
# ----------------------------------------------------------------------


function build_selected_surfaces(state, input, selected_steps)
    comp = state.compaction
    nx, ny = size(comp.n_layers)

    surfaces = Vector{Matrix{Float32}}(undef, length(selected_steps))

    for (idx, t) in enumerate(selected_steps)
        surf = zeros(Float32, nx, ny)

        for i in 1:nx, j in 1:ny
            z = 0.0f0

            for k in reverse(1:comp.n_layers[i,j])
                step_k = comp.layer_step[i,j,k]
                step_k > t && continue

                fidx = comp.layer_facies[i,j,k]
                H0   = comp.layer_thickness0[i,j,k]

                facies = input.facies[fidx]
                ϕ0 = facies.depositional_porosity

                z_mid = z + 0.5f0 * H0
                ϕ = Production.porosity_at_depth(facies, z_mid)
                Hc = Production.compacted_thickness(H0, ϕ0, ϕ)

                z += Hc
            end

            surf[i,j] = z
        end

        surfaces[idx] = surf
    end

    return surfaces
end

end