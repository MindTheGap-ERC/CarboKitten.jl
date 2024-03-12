# Unitful
Physical quantities in CarboKitten are always specified using the `Unitful.jl` framework.

``` {.julia file=test/Unitful.jl}
@testset "Unitful" begin
    using Unitful
    using Unitful.DefaultSymbols

    <<unitful-spec>>
end
```

## Variables vs. string macros
`Unitful` package offers two basic ways to enter quantities: either using predefined symbols (polluting your namespace with one-letter variables), or using special string macros.

``` {.julia #unitful-spec}
@test 1.0m === 1.0u"m"
@test 42J/s == 42u"W"
```

In many cases, your CarboKitten scripts will contain little else than the input specification. In such a case `using Unitful.DefaultSymbols` gives a bit cleaner, more readable look.

## Reading specs
Suppose we simulate a pendulum. We would have an input spec defined as follows:

``` {.julia #unitful-spec}
@kwdef struct Pendulum
    length :: typeof(1.0m)
    time_step :: typeof(1.0s)
    phi0 :: typeof(1.0rad)
    omega0 :: typeof(1.0rad/s)
end
```

Then input can be given as follows:

``` {.julia #unitful-spec}
pendulum = Pendulum(
    length = 2.0m,
    time_step = 1ms,
    phi0 = 30°,
    omega0 = 0rad/s
)
```

Units are automatically converted to the types specified in the API.

``` {.julia #unitful-spec}
@test pendulum.time_step === 0.001s
@test pendulum.phi0 === (π/6)rad
```

## Dimensions
Unitful has dimensions of length `𝐋`, mass `𝐌` and time `𝐓` as bold upper-case Unicode symbols. These can be entered in Julia with `\bfL`, `\bfM` etc.
When you define a function that needs, say, an energy, which has SI units of ${\rm J = (m/s)^2\ kg}$, we can construct the dimensions. Defining a few constants:

``` {.julia #unitful-spec}
let 𝐄 = (𝐋/𝐓)^2 * 𝐌,
    h = 6.62607015e-34u"J*s",
    c = 299792458u"m/s"
    <<unitful-photon-example>>
end
```

We can abstract over the specific units by defining a generic method. Now we can compute the wavelength of a photon, given its energy in any unit of energy.

``` {.julia #unitful-photon-example}
photon_wave_length(E::Quantity{Float64,𝐄,J}) where {J} =
    uconvert(u"Å", h * c / E)

@test photon_wave_length(2.38u"eV") ≈ 5209.4201u"Å"
@test_throws MethodError photon_wave_length(1u"m")
```
