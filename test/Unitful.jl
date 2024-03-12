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
    # ~/~ begin <<docs/src/unitful.md#unitful-spec>>[4]
    let 𝐄 = (𝐋/𝐓)^2 * 𝐌,
        h = 6.62607015e-34u"J*s",
        c = 299792458u"m/s"
        # ~/~ begin <<docs/src/unitful.md#unitful-photon-example>>[init]
        photon_wave_length(E::Quantity{Float64,𝐄,J}) where {J} =
            uconvert(u"Å", h * c / E)

        @test photon_wave_length(2.38u"eV") ≈ 5209.4201u"Å"
        @test_throws MethodError photon_wave_length(1u"m")
        # ~/~ end
    end
    # ~/~ end
end
# ~/~ end