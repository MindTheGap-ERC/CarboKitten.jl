# ~/~ begin <<docs/src/memory-writer.md#examples/autocycles.jl>>[init]
module AutoCycles

using CarboKitten
using CarboKitten.Models: WithoutCA as M
using CarboKitten.Visualization: profile_plot!

using Makie
using CarboKitten: Box

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

function run()
    CarboKitten.init()

    facies = [
        M.Facies(
            maximum_growth_rate=100.0u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=100.0u"m/yr",
            wave_velocity=v_const(-5.0u"m/yr")),
        M.Facies(
            maximum_growth_rate=20.0u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=10.0u"m/yr",
            wave_velocity=v_const(0.0u"m/yr"))]

    
    input = M.Input(
        box=CarboKitten.Box{Coast}(grid_size=(500, 1), phys_scale=50.0u"m"),
        time=TimeProperties(
            Δt=20u"yr",
            steps=40000),

        output=Dict(
            :profile => OutputSpec(slice = (:, 1), write_interval = 40),
            :col1 => OutputSpec(slice = (100, 1)),
            :col2 => OutputSpec(slice = (200, 1)),
            :col3 => OutputSpec(slice = (300, 1)),
            :col4 => OutputSpec(slice = (400, 1))),

        initial_topography=(x, y) -> (15u"km" - x) / 300.0,
        sea_level=t -> 0.0u"m",

        subsidence_rate=50.0u"m/Myr",
        disintegration_rate=100.0u"m/Myr",
        precipitation_time=50.0u"yr",

        insolation=400.0u"W/m^2",
        sediment_buffer_size=50,
        depositional_resolution=0.5u"m",
        facies=facies,

        transport_solver=Val{:forward_euler},
        intertidal_zone=0.0u"m")

    result = run_model(Model{M}, input, MemoryOutput(input))
end

function plot(result::MemoryOutput)
    header = result.header
    slice = result.data_slices[:profile]
    n_cols = length(result.data_columns)

	fig = Figure()
    ax1 = Axis(fig[1, 1:n_cols])

	x = header.axes.x
	t = header.axes.t
	
	plot = profile_plot!(ax1, header, slice; colorrange=(0.2, 1.0)) do x; x[1] / sum(x) end 
    col_positions = [x[col.slice[1]] |> in_units_of(u"km") for col in values(result.data_columns)]
    vlines!(ax1, col_positions; color=:red)

    Colorbar(fig[1, n_cols+1], plot; label=L"f_1 / f_{\textrm{total}}")

    col_names = sort!(collect(keys(result.data_columns)))
    for (i, k) in enumerate(col_names)
        f1 = result.data_columns[k].deposition[1,:]
        f2 = result.data_columns[k].deposition[2,:]
        f_total = f1 .+ f2
        ax = Axis(fig[2, i], title=string(k) * " amount")
        lines!(ax, f1 |> in_units_of(u"m"), t[1:end-1] |> in_units_of(u"Myr"); label="facies 1")
        lines!(ax, f2 |> in_units_of(u"m"), t[1:end-1] |> in_units_of(u"Myr"); label="facies 2")
        lines!(ax, f_total |> in_units_of(u"m"), t[1:end-1] |> in_units_of(u"Myr"); label="total", color=:black, linewidth=2)
        axislegend(ax)
    end

	fig
end

end

# ~/~ end
