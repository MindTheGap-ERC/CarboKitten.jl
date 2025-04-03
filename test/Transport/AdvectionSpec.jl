# ~/~ begin <<docs/src/finite-difference-transport.md#test/Transport/AdvectionSpec.jl>>[init]
@testset "CarboKitten.Transport.Advection.scale-invariance" begin

using CarboKitten
using CarboKitten.Transport.Solvers: runge_kutta_4
using CarboKitten.Transport.Advection: transport
using Unitful

let box = Box{Periodic{2}}(grid_size=(32, 32), phys_scale=1.0u"m")
    solver = runge_kutta_4(Float64, box)
    wave_velocity = _ -> ((0.5u"m/s", 0.0u"m/s"), (0.0u"1/s", 0.0u"1/s"))
    diffusivity = 5.0u"m/s"
    w = randn(box.grid_size...) * u"m"
    C1 = randn(box.grid_size...)
    C2 = C1 .* 10.0
    dt = 1.0u"s"
    df(C, _) = transport(box, diffusivity, wave_velocity, C, w)

    for i = 1:10
        solver(df, C1, 0.0u"s", dt)
        solver(df, C2, 0.0u"s", dt)
    end

    @test isapprox(C1 .* 10.0, C2)
end

end
# ~/~ end
