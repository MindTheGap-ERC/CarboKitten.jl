# Carbonate Production

The paper by Bosscher and Schlager (1992) [Bosscher1992](@cite) is an early computer model for simulating reef growth. This paper contains some of the essential ingredients that we find back in CarboKitten. We reproduce their results within the framework of CarboKitten's larger design.

The BS92 model assumes a direct relation between water depth and sediment accumulation rate. That way we can model reef growth by integrating an Ordinary Differential Equation (ODE). The [`Production` component](./components/production.md) provides this model for the rest of CarboKitten.

![summary plot](fig/bs92-summary.png)

## Parameters

> * Maximum growth rate $G_m$. The maximum rate of reef growth is in the range of $10-15 {\rm mm\ yr^{-1}}$ (Macintyre etal., 1977; Adey, 1978; Davies, 1983).
> * Extinction coefficient $k$. This is a measure of the extinction of photosynthetically active radiation (PAR), i.e. light with a wavelength of 400-700 nm.  The value of k for oceanic waters ranges from $0.04$ to $0.16\ {\rm m^{-1}}$ (Jerlov, 1976); reported values for reef waters also lie within this range (Brakel, 1979; Van den Hoek et al., 1975; Weinberg, 1976; Chalker, 1981; Porter, 1985).
> * Surface light intensity $I_0$. The light intensity at the water surface at midday in the tropics lies in the range of $2000-2250\ {\rm \mu E\ m^{-2}s^{-1}}$.
> * Saturating light intensity $I_k$. Light saturating intensities are in the range $50-450\ {\rm \mu E\ m^{-2}s^{-1}}$, depending on species and water depth (Chalker, 1981; Wyman et al., 1987). Photoadaptation of reef-building corals has not been taken into account. More generally, light does not become a limiting factor for coral growth until it reaches roughly 10% of its surface value (B. E. Chalker, in Done, 1983).
*from BS92*

## Growth Rate

The growth rate is given as

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right),$$

``` {.julia #carbonate-production}
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))
```

``` {.julia #b92-model}
<<carbonate-production>>

struct Parameters
     I₀::Float64
     Iₖ::Float64
     k::Float64
     gₘ::Float64
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)
```

where $w$ is the water depth in meters, $g_m$ is the maximum growth rate in ${\rm m}\ {\rm My}^-1$, $I_0$ is surface light intensity, $I_k$ is saturation light intensity, and $k$ is the extinction coefficient. We have exponential decay of light intensity as we get to deeper water, and the carbonate factories respond to light intensity through a $\tanh$ (hyperbolic tangent function). This is by no means an exact relation, rather the $\tanh$ function interpolates smoothly between one and zero. We specify a maximum growth rate and a typical intensity at which the species is no longer productive.

The shape of $\tanh \circ \exp$ look like this:

![Tangens hyperbolicus](fig/tanh.svg)

Notice that the numbers inside the exponential need to be unit-free, so does the output. The value of $\tanh \circ \exp$ at a depth of 0 is $0.7615\dots$. This does not make much sense, as I believe we should start at a value of 1 at the surface. By setting $I_0 / I_k$ to some value $>1$ this can be alleviated, but it changes the interpretation of the constants a little. The idea is that above a certain insolation, light is not the limiting factor to the rate of photosynthesis.

To reproduce Figure 2 in B13, I had to change the values for $g_m$ to 500, 250, and 125 respectively, the other values from Table 2 remained the same. I guess this was done for illustration purposes.

![Production curves](fig/b13-fig2.svg)

```@raw html
<details><summary>Plotting code</summary>
```

``` {.gnuplot .build file=examples/plot-tanh.gnuplot target=docs/src/_fig/tanh.svg}
set term svg size 700, 300 font "sans serif, 14" linewidth 1.5
set xrange [-5:10]
set yrange [-0.1:1.1]
set grid
set key outside
set xlabel "x"
set ylabel "y"
plot tanh(exp(-x)) lc rgb 'black', tanh(exp(4)*exp(-x)), tanh(exp(-0.5*x))
```

``` {.gnuplot .build file=examples/burgess2013-fig2.gnuplot target=docs/src/_fig/b13-fig2.svg}
set term svg size 500, 600 font "sans serif,14" linewidth 1.5
set trange [0:100]
set yrange [100:0]
set xrange [-20:520]
set parametric
set key right bottom
set grid
set ylabel "Water depth (m)"
set xlabel "Production rates"
plot 500*tanh(6.7 * exp(-0.8 * t)), t title 'Carbonate factory 1', \
     250*tanh(6.7 * exp(-0.1 * t)), t title 'Carbonate factory 2', \
     125*tanh(6.7 * exp(-0.005 * t)), t title 'Carbonate factory 3'
```

```@raw html
</details>
```

## Depth Evolution

The use of water depth in both BS92 and B13 can be a bit confusing. Plots are shown up-side-down and little is done to disambiguate depth with sea level rising or lowering, or sediment accreting. Growth in deposition should give shallower sea bed. BS92 write $w = (h_0 + h(t)) - (s_0 + s(t))$. Actually $s_0$ is best set to 0, or simply included into $s(t)$ and $h_0$ can be replaced with setting $h(t=0) = h_0$. Then, as we have the growth rate as a function of water depth $g(w)$, we can say

$$\partial_t h = -g_m {\rm tanh}\left[\frac{I_0}{I_k} \exp(-k (h - s(t)))\right].$$

``` {.julia #b92-model}
function model(p::Parameters, s, t_end::Float64, h₀::Float64)
     ∂h(h::Float64, t::Float64) = let w = h - s(t)
          w >= 0.0 ? -g(p, h - s(t)) : 0.0
     end

     dt = 1000.0
     times = 0.0:dt:t_end
     result = zeros(Float64, length(times))
     result[1] = h₀
     for (i, t) in enumerate(times[1:end-1])
          h = result[i]
          for j = 0:99
               h += ∂h(h, t + j * dt/100) * (dt/100)
          end
          result[i+1] = h
     end
     return result
end
```

It seems Eq. 5 in BS92 (the most important equation in the paper mind you!) is missing both a minus sign and a set of parentheses. Also, we should remark that at negative depth (subareal exposure) we should halt all growth.

## Crosssection

The most impressive result in BS92 is the last figure. They show an input curve for $s(t)$ but give no functional description. The curve starts with a linear drop from 0 to 120m depth over a time of 20000 years, then slowly rises with $s(t) = a +  bt + A \sin(2\pi t / P)$, with a period $P = \sim 15-20 {\rm kyr}$, amplitude $A = \sim 40 {\rm m}$. It might be easiest to take a screenshot of the PDF and convert the graph into a table.

Using `DifferentialEquations.jl` we can integrate Equation @eq:growth-eqn. Interestingly, the only integrator that gave me noise free results is `Euler`. This may be due to the sudden shut-down of production at $w = 0$.

``` {.julia file=examples/model/bs92/using_ode.jl}
module BS92

using CarboKitten.DataSets: bosscher_schlager_1992
using Interpolations
using Unitful

<<b92-model>>

function sealevel_curve()
     data = bosscher_schlager_1992()
     linear_interpolation(data.time / u"yr", data.sealevel / u"m")
end

struct Scenario
     param::Parameters
     sealevel
     t_end::Float64
end

model(s::Scenario, h₀::Float64) = model(s.param, s.sealevel, s.t_end, h₀)

SCENARIO_A = Scenario(
     Parameters(2000.0, 250.0, 0.05, 0.005),
     sealevel_curve(),
     80_000.0)

end
```

Finally, we can try to reproduce figure 8 in BS92.

![stratigraphy following fig. 8 in BS92](fig/bs92-fig8.svg)

Note the simplicity of this result: there is no dependency on space, only on the initial depth $h_0$.

```@raw html
<details><summary>Plotting code</summary>
```

``` {.julia .task file=examples/model/bs92/fig8.jl}
#| creates: docs/src/_fig/bs92-fig8.svg
#| requires: examples/model/bs92/using_ode.jl
#| collect: figures

module Script
     include("using_ode.jl")
     using CairoMakie

     function main()
          h0 = LinRange(0, 200, 101)
          result = hcat([BS92.model(BS92.SCENARIO_A, h) for h in h0]...)
          t = LinRange(0, 80_000, 81)

          fig = Figure(resolution=(600,900))
          ax = Axis(fig[1,1], xlabel="initial depth (m)", ylabel="depth (m)", yreversed=true)
          for l in eachrow(result)
               lines!(ax, h0, vec(l); color=:steelblue4, linewidth=0.5)
          end
          ax = Axis(fig[2,1], xlabel="time (years)", ylabel="depth (m)", yreversed=true)
          lines!(ax, t, BS92.SCENARIO_A.sealevel(t); color=:steelblue4)

          save("docs/src/_fig/bs92-fig8.svg", fig)
     end
end

Script.main()
```

```@raw html
</details>
```

## BS92 in CarboKitten stack

```component-dag
CarboKitten.Models.BS92
```

Within the CarboKitten design, we can express the BS92 model a bit more succinctly. The following produces output that is fully compatible with other CarboKitten models and the included post processing and visualization stack. The `H5Writer` module provides a `run` method that expects the `initial_state`, `step!` and `write_header` methods to be available.

``` {.julia file=src/Models/BS92.jl}
@compose module BS92
@mixin Tag, H5Writer, Production

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    return State(0, sediment_height)
end

function step!(input::Input)
    τ = uniform_production(input)
    function (state::State)
        prod = τ(state)
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
        state.step += 1
        return Frame(
            production = prod,
            deposition = prod)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

end
```

``` {.julia .task file=examples/model/bs92/run.jl}
#| creates: data/output/bs92.h5

module Script

using CarboKitten
using CarboKitten.DataSets: bosscher_schlager_1992

using Interpolations
using Unitful

function sealevel_curve()
    data = bosscher_schlager_1992()
    linear_interpolation(data.time, data.sealevel)
end

const INPUT = BS92.Input(
    tag = "example model BS92",
    box = Box{Coast}(grid_size=(100, 1), phys_scale=600.0u"m"),
    time = TimeProperties(
      Δt = 10.0u"yr",
      steps = 8000,
      write_interval = 100),
    sea_level = let sc = sealevel_curve()
      t -> -sc(t)
    end,
    initial_topography = (x, y) -> - x / 300.0,
    subsidence_rate = 0.0u"m/yr",
    insolation = 400.0u"W/m^2",
    facies = [BS92.Facies(
      maximum_growth_rate = 0.005u"m/yr",
      saturation_intensity = 50.0u"W/m^2",
      extinction_coefficient = 0.05u"m^-1"
    )])

function main()
    run_model(Model{BS92}, INPUT, "data/output/bs92.h5")
end

end

Script.main()
```

``` {.julia .task file=examples/model/bs92/plot.jl}
#| creates: docs/src/_fig/bs92-summary.png
#| requires: data/output/bs92.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/bs92-summary.png", summary_plot("data/output/bs92.h5"))
```

## Multiple facies

Using the above implementation of the model by Bosscher and Schlager, we can run the same model with multiple facies. We use the same parameters as Burgess2013, but divide the `maximum_growth_rate` by four, since we have no CA running and all facies produce at the same time.


``` {.julia .task file=examples/model/bs92/multi-facies-run.jl}
#| creates: data/output/bs92-multi-facies.h5

module Script

using CarboKitten

const FACIES = [
    BS92.Facies(
         maximum_growth_rate=500u"m/Myr"/4,
         extinction_coefficient=0.8u"m^-1",
         saturation_intensity=60u"W/m^2"),
    BS92.Facies(
         maximum_growth_rate=400u"m/Myr"/4,
         extinction_coefficient=0.1u"m^-1",
         saturation_intensity=60u"W/m^2"),
    BS92.Facies(
         maximum_growth_rate=100u"m/Myr"/4,
         extinction_coefficient=0.005u"m^-1",
         saturation_intensity=60u"W/m^2")]

const INPUT = BS92.Input(
    tag = "example model BS92",
    box = Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
    time = TimeProperties(
        Δt = 200.0u"yr",
        steps = 5000,
        write_interval = 1),
    sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    initial_topography = (x, y) -> - x / 300.0,
    subsidence_rate = 50.0u"m/Myr",
    insolation = 400.0u"W/m^2",
    facies = FACIES)

function main()
    run_model(Model{BS92}, INPUT, "data/output/bs92-multi-facies.h5")
end

end

Script.main()
```

``` {.julia .task file=examples/model/bs92/multi-facies-plot.jl}
#| creates: docs/src/_fig/bs92-multi-facies.png
#| requires: data/output/bs92-multi-facies.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/bs92-multi-facies.png", summary_plot("data/output/bs92-multi-facies.h5", wheeler_smooth=(3, 5)))
```

![BS92 with multiple facies](fig/bs92-multi-facies.png)
