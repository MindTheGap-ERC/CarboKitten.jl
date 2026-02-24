# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/disintegration-transfer.jl>>[init]

module DisintTransferExample

using CarboKitten
using CarboKitten: Box
using CarboKitten.Models: WithoutCA as M
using CarboKitten.Visualization: profile_plot!
using CairoMakie

# for setting a constant wave velocity, if required
v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

# ~/~ begin <<docs/src/active-layer-transport.md#disint-transfer-example>>[init]
# put in some sediment to transport
function initial_sediment(x,y)
    if x > 10u"km" && x < 15u"km"
        return 5.0u"m"
    else
        return 0.0u"m"
    end
end

function run()
    facies = [
        M.Facies(
            maximum_growth_rate=0.0u"m/Myr",
            initial_sediment=initial_sediment),
        M.Facies(
            diffusion_coefficient=100.0u"m/yr",
            wave_velocity=v_const(0.0u"m/yr")),
    ]

    input = M.Input(
        box=Box{Coast}(grid_size=(100, 1), phys_scale=250.0u"m"),
        time=TimeProperties(
            Δt=20.0u"yr",
            steps=40000),
        output=Dict(:full=>OutputSpec(write_interval=10),
                    :profile => OutputSpec(slice=(:,1), write_interval=40)),
        initial_topography=(x, y) -> -25u"m",
        sea_level=t -> 0.0u"m",
        subsidence_rate=0.0u"m/Myr",
        disintegration_rate=10.0*u"m/Myr",
        cementation_time=10.0*u"yr",
        insolation=400.0u"W/m^2",
        sediment_buffer_size=50,
        depositional_resolution=0.5u"m",
        facies=facies,
        transport_solver=Val{:forward_euler},
        intertidal_zone=0.0u"m",
        disintegration_transfer=f->stack((0.0.*f[1,:,:], f[1,:,:].+f[2,:,:]), dims=1),
        save_active_layer=true)

    result = run_model(Model{M}, input, MemoryOutput(input))
    return result
end
# ~/~ end

function main()
    res = run()
    data = res.data_slices[:profile]
    n_facies = res.header.n_facies
    fig = Figure()
    ax = Axis(fig[1,1])
    profile_plot!(argmax, ax, res.header, data; alpha=1.0, colormap=cgrad(Makie.wong_colors()[1:n_facies], n_facies, categorical=true))
    save("docs/src/_fig/disintegration-transfer-example.png",fig)
end

end
DisintTransferExample.main()
# ~/~ end
