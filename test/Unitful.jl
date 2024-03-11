# ~/~ begin <<docs/src/unitful.md#test/Unitful.jl>>[init]
@testset "Unitful" begin
    using Unitful
    using Unitful.DefaultSymbols

    # ~/~ begin <<docs/src/unitful.md#unitful-spec>>[init]
    @test 1.0m === 1.0u"m"
    @test 42J/s == 42u"W"
    # ~/~ end
    # ~/~ begin <<docs/src/unitful.md#unitful-spec>>[1]
    @kwdef struct Pendulum
        length :: typeof(1.0m)
        time_step :: typeof(1.0s)
        phi0 :: typeof(1.0rad)
        omega0 :: typeof(1.0rad/s)
    end
    # ~/~ end
    # ~/~ begin <<docs/src/unitful.md#unitful-spec>>[2]
    pendulum = Pendulum(
        length = 2.0m,
        time_step = 1ms,
        phi0 = 30°,
        omega0 = 0rad/s
    )
    # ~/~ end
    # ~/~ begin <<docs/src/unitful.md#unitful-spec>>[3]
    @test pendulum.time_step === 0.001s
    @test pendulum.phi0 === (π/6)rad
    # ~/~ end
end
# ~/~ end