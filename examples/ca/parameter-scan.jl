# ~/~ begin <<docs/src/carbocat-ca.md#examples/ca/parameter-scan.jl>>[init]
module Script

using CarboKitten
using CarboKitten.Burgess2013
using CarboKitten.Stencil
using CarboKitten.Utility
using GLMakie
using .Iterators: peel, drop

function main()
    fig = Figure(resolution=(2000, 2000))
    for i in 4:12
        for j in (i+1):12
            print(".")
            gl = fig[i, j] = GridLayout()
            for k in i:j
                # for l in k:j
                let l = j
                    init = rand(0:3, 50, 50)
                    facies = [
                        Facies((i, j), (k, l), 0, 0, 0),
                        Facies((i, j), (k, l), 0, 0, 0),
                        Facies((i, j), (k, l), 0, 0, 0),
                    ]
                    (result, _) = peel(drop(run_ca(Periodic{2}, facies, init, 3), 10))

                    ax = Axis(gl[k-i, l-k], aspect=AxisAspect(1),
                              xticksvisible=false, xticklabelsvisible=false,
                              yticksvisible=false, yticklabelsvisible=false)
                    heatmap!(ax, result)
                end
            end
        end
    end
    save("docs/src/fig/parameter-scan.png", fig)
end

end  # module Script

Script.main()
# ~/~ end