module CAExplorer

using CarboKitten
using CarboKitten.Burgess2013
using CarboKitten.Stencil
using CarboKitten.Utility
using GLMakie
using Observables
using .Iterators: peel, drop

function tic(t, dt, running)
    @async begin
        while running[]
            sleep(dt)
            t[] += dt
        end
    end
end

function main()
    fig = Figure(resolution=(800, 800))
    init = rand(0:3, 50, 50)
    image = Observable(init)

    ax = Axis(fig[1, 1], aspect=AxisAspect(1),
              xticksvisible=false, xticklabelsvisible=false,
              yticksvisible=false, yticklabelsvisible=false)

    sg = SliderGrid(fig[2, 1],
        ( label = "min viability", range = 0:25, startvalue = 4 ),
        ( label = "max viability", range = 0:25, startvalue = 10 ),
        ( label = "min activation", range = 0:25, startvalue = 6 ),
        ( label = "max activation", range = 0:25, startvalue = 10 ))

    # add play button
    t = Observable(0.0)
    fig[3, 1] = button_grid = GridLayout(tellwidth = false)
    reset_button = button_grid[1, 1] = Button(fig, label="Reset")
    step_button = button_grid[1, 2] = Button(fig, label="Step")
    running = button_grid[1, 3] = Toggle(fig, active=false)
    button_grid[1, 4] = Label(fig, lift(x -> x ? "Playing" : "Paused", running.active))

    on(running.active) do activated
        if activated
            tic(t, 0.05, running.active)
        end
    end

    ca = lift([s.value for s in sg.sliders]..., reset_button.clicks) do i, j, k, l, _
        facies = [
            Facies((i, j), (k, l), 0, 0, 0),
            Facies((i, j), (k, l), 0, 0, 0),
            Facies((i, j), (k, l), 0, 0, 0),
        ]
        run_ca(Periodic{2}, facies, copy(init), 3)
    end

    onany(t, step_button.clicks) do _, _
        image[] = take!(ca[])
    end

    heatmap!(ax, image)
    fig
end

end # module CAExplorer