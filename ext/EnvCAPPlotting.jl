# ~/~ begin <<docs/src/visualization/envcap_plotting.md#ext/EnvCAPPlotting.jl>>[init]
module EnvCAPPlotting

using CarboKitten
using CarboKitten.Models.EnvCAP.EnvMapping: dominant_env_block, env_to_factory_prior_block
using CarboKitten.Export: read_volume
using CairoMakie
using Unitful

export plot_envcap_maps_at_z, plot_envcap_run_outputs

const DEFAULT_ENV_TO_FACTORY = [
    0.8  0.15  0.05;
    0.3  0.55  0.15;
    0.05 0.25  0.70
]

const DEFAULT_TAGS = [
    "envcap_ref0",
    "envcap_ref05",
    "envcap_ref1",
]

function as_meters(a)
    try
        return Float64.(ustrip.(u"m", a))
    catch
        try
            return Float64.(ustrip.(a))
        catch
            return Float64.(a)
        end
    end
end

function checked_z_index(block, z_index; name)
    nz = size(block, ndims(block))

    if z_index > nz
        @warn "$name has only $nz z layers. Plotting z=$nz instead of z=$z_index."
        return nz
    elseif z_index < 1
        @warn "$name received z=$z_index. Plotting z=1 instead."
        return 1
    end

    return z_index
end

function factory_block_from_deposition(deposition; dz = 0.5, min_nz = 1)
    dep = as_meters(deposition)

    n_factories, nx, ny, nt = size(dep)

    total_thickness = dropdims(sum(dep; dims = (1, 4)), dims = (1, 4))
    nz = max(min_nz, 1, ceil(Int, maximum(total_thickness) / dz))

    factory_block = zeros(Int, nx, ny, nz)

    for ix in 1:nx, iy in 1:ny
        zpos = 1

        for it in 1:nt
            column = dep[:, ix, iy, it]
            thickness = sum(column)

            thickness <= 0.0 && continue

            factory = argmax(column)
            nvox = max(1, round(Int, thickness / dz))
            ztop = min(nz, zpos + nvox - 1)

            factory_block[ix, iy, zpos:ztop] .= factory

            zpos = ztop + 1
            zpos > nz && break
        end
    end

    return factory_block
end

function read_shared_conditioning(
        stage1_file;
        env_to_factory = DEFAULT_ENV_TO_FACTORY,
        dz = 0.5,
        min_nz = 1)

    _, stage1_data = read_volume(stage1_file, :full)

    environment_belt = dominant_env_block(
        stage1_data.deposition;
        dz = dz,
        min_nz = min_nz,
    )

    factory_prior = env_to_factory_prior_block(
        environment_belt,
        env_to_factory,
    )

    return environment_belt, factory_prior
end

function add_map_panel!(
        parent,
        row,
        col,
        data;
        title,
        colorrange,
        colorbar_label,
        colormap = :viridis,
        categorical = false,
        category_ticks = nothing,
        map_width = 320,
        map_height = 220)

    panel = GridLayout()
    parent[row, col] = panel

    ax = Axis(
        panel[1, 1],
        title = title,
        xlabel = "x index",
        ylabel = "y index",
        aspect = DataAspect(),
        titlesize = 13,
        xlabelsize = 11,
        ylabelsize = 11,
        xticklabelsize = 9,
        yticklabelsize = 9,
    )

    hm = heatmap!(
        ax,
        data;
        colorrange = colorrange,
        colormap = colormap,
    )

    if categorical && category_ticks !== nothing
        Colorbar(
            panel[1, 2],
            hm;
            label = colorbar_label,
            width = 14,
            labelsize = 10,
            ticklabelsize = 9,
            ticks = category_ticks,
        )
    else
        Colorbar(
            panel[1, 2],
            hm;
            label = colorbar_label,
            width = 14,
            labelsize = 10,
            ticklabelsize = 9,
        )
    end

    colgap!(panel, 8)
    colsize!(panel, 1, Fixed(map_width))
    colsize!(panel, 2, Fixed(50))
    rowsize!(panel, 1, Fixed(map_height))

    return ax
end

function plot_envcap_maps_at_z(
        tag;
        output_dir = "data/output",
        stage1_file = joinpath(output_dir, "envcap_stage1.h5"),
        envcap_file = joinpath(output_dir, "$(tag).h5"),
        env_to_factory = DEFAULT_ENV_TO_FACTORY,
        z_index = 41,
        dz = 0.5,
        fig_dir = joinpath(output_dir, "fig"))

    mkpath(fig_dir)

    environment_belt, factory_prior = read_shared_conditioning(
        stage1_file;
        env_to_factory = env_to_factory,
        dz = dz,
        min_nz = z_index,
    )

    _, envcap_data = read_volume(envcap_file, :full)

    factory_block = factory_block_from_deposition(
        envcap_data.deposition;
        dz = dz,
        min_nz = z_index,
    )

    z_belt = checked_z_index(environment_belt, z_index; name = "environment_belt")
    z_prior = checked_z_index(factory_prior, z_index; name = "factory_prior")
    z_factory = checked_z_index(factory_block, z_index; name = "factory_block")

    belt_map = environment_belt[:, :, z_belt]
    factory_map = factory_block[:, :, z_factory]

    n_envs = size(env_to_factory, 1)
    n_factories = size(factory_prior, 1)
    z_m = round((z_index - 1) * dz; digits = 2)
    
    # Categorical colormaps
    env_cmap = cgrad(:tab10, n_envs; categorical = true)
    fac_cmap = cgrad(:tab10, n_factories; categorical = true)
    
    fig = Figure(
        size = (1900, 950),
        fontsize = 12,
        figure_padding = (35, 35, 35, 35),
    )
    
    Label(
        fig[1, 1:2],
        "$(tag): environment belt, probability prior, and deposited factory at z=$(z_index) ≈ $(z_m) m",
        fontsize = 18,
        halign = :left,
        tellwidth = false,
    )
    
    left = GridLayout()
    right = GridLayout()
    
    fig[2, 1] = left
    fig[2, 2] = right
    
    colgap!(fig.layout, 35)
    rowgap!(fig.layout, 20)
    colgap!(left, 20)
    rowgap!(left, 20)
    
    # ---- Left 2x2 block ----
    
    add_map_panel!(
        left,
        1, 1,
        belt_map;
        title = "Environment belt\nz=$(z_belt)",
        colorrange = (0.5, n_envs + 0.5),
        colorbar_label = "Environment",
        colormap = env_cmap,
        categorical = true,
        category_ticks = 1:n_envs,
        map_width = 360,
        map_height = 240,
    )
    
    if n_factories >= 1
        add_map_panel!(
            left,
            1, 2,
            factory_prior[1, :, :, z_prior];
            title = "P(factory 1)\nz=$(z_prior)",
            colorrange = (0.0, 1.0),
            colorbar_label = "Probability",
            colormap = :viridis,
            map_width = 360,
            map_height = 240,
        )
    end
    
    if n_factories >= 2
        add_map_panel!(
            left,
            2, 1,
            factory_prior[2, :, :, z_prior];
            title = "P(factory 2)\nz=$(z_prior)",
            colorrange = (0.0, 1.0),
            colorbar_label = "Probability",
            colormap = :viridis,
            map_width = 360,
            map_height = 240,
        )
    end
    
    if n_factories >= 3
        add_map_panel!(
            left,
            2, 2,
            factory_prior[3, :, :, z_prior];
            title = "P(factory 3)\nz=$(z_prior)",
            colorrange = (0.0, 1.0),
            colorbar_label = "Probability",
            colormap = :viridis,
            map_width = 360,
            map_height = 240,
        )
    end
    
    # ---- Right large factory map ----
    
    add_map_panel!(
        right,
        1, 1,
        factory_map;
        title = "Deposited factory\nz=$(z_factory)",
        colorrange = (0.5, n_factories + 0.5),
        colorbar_label = "Factory",
        colormap = fac_cmap,
        categorical = true,
        category_ticks = 1:n_factories,
        map_width = 520,
        map_height = 520,
    )
    
    resize_to_layout!(fig)
    
    outfile = joinpath(fig_dir, "$(tag)_z$(z_index)_maps.png")
    save(outfile, fig)
    
    println("Saved $(outfile)")
    return fig
end

function plot_envcap_run_outputs(;
        tags = DEFAULT_TAGS,
        output_dir = "data/output",
        env_to_factory = DEFAULT_ENV_TO_FACTORY,
        z_index = 41,
        dz = 0.5)

    for tag in tags
        plot_envcap_maps_at_z(
            tag;
            output_dir = output_dir,
            env_to_factory = env_to_factory,
            z_index = z_index,
            dz = dz,
        )
    end

    return nothing
end

function main()
    plot_envcap_run_outputs()
    return nothing
end

end
# ~/~ end
