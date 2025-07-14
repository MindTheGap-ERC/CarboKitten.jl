# Active Layer Transport

The following is inspired on well-known **active layer** approaches in river bed sediment transport [Paola1992](@cite) [James2010](@cite) [^1]. All quantities with subscript $f$ are facies dependent. Sediment is measured in meters of deposited material. $P_f$ is the production of sediment per facies in $m/s$. Further unit calculations would be more readable if we consider the unit of sediment as separate, so for instance it doesn't cancel against $m^2$ in the units of sediment flux. In the implementation, $\nu$ has the units of ${\rm m}$ which is totaly weird. TBC

[^1]: Literature on active (or mixing) layer transport modeling is vast. Most of which is concerned with much smaller time scales, and more complicated physics than we are mostly dealing with.

In a model without transport, we could write

$$\sigma + \sum_f {{\partial \eta_f} \over {\partial t}} = \sum_f P_f,$$

where $\sigma$ is the subsidence rate in $m/s$. We consider the mass balance for each facies separately.

We suppose that loose sediment, either fresh production or disintegrated older sediment, is being transported in a layer on top of the sea bed. The flux in this layer is assumed to be directly proportional to the local slope of the sea bed $| \nabla_x \eta_* |$, where $\eta_* = \sum_f \eta_f$, the sum over all facies contributions, including $\eta_0$, the initial bedrock eleveation.

![Schematic of Active Layer approach](fig/active-layer-export.svg)

The active layer now contains a concentration $C_f$ particles of different grain size (for each facies $f$). If needed, $C_f = \alpha_f P_f$ where $\alpha_f$ is some facies parameter determining the fraction of production that is available for transport. The sediment flux is given as,

$${\bf q_f} = -\nu_f C_f {\bf \nabla_x} \eta_*.$$

The following is the mass balance:

$$\sigma + {{\partial \eta_*} \over {\partial t}} = -\sum_f {\bf \nabla_x} \cdot {\bf q_f} + \sum_f P_f,$$

In our modelling we keep track of individual contributions per facies over time [^2].

[^2]: Note that in other approaches to active layer transport, like Paola 1992, there would be a factor $1/C_f$. Here we have a different interpretation to what the concentration means: the sediment settles down after transport, such that the concentration has no impact on the change in sediment surface elevation.

Combining these equations, and ignoring subsidence for the moment (which is a global effect and can't be expressed on a per-facies basis), we get a component-wise diffusion equation

$${{\partial \eta_f(x)}\over{\partial t}} = {\bf \nabla_x} \cdot \big[ \nu_f \alpha_f\ P_f(x)\ {\bf \nabla_x} \eta_{*}(x) \big] + P_f(x),$$

In our model we need to solve this equation one time-step each iteration. If we solve this using forward methods, we should be reminded of the CFL limit for diffusion equations (depending on the diffusion constants and grid size we shouldn't pick the time steps too large). Alternatively, for these two-dimensional situations, an implicit approach is feasible. Also we should take care that somehow $\nabla(\nu\alpha P \nabla \eta) + P > 0$. The interpretation being that we can't transport more than we produce, even if there is capacity to do so.

To solve this equation, it is nicer to expand the transport-diffusion term using the product rule, in short notation:

$$\partial_t \eta_f = \nu' \nabla P_f(x) \cdot \nabla \eta(x) + \nu' P_f(x) \nabla^2 \eta(x) + P_f(x),$$

where $\nu' = \nu_f \alpha_f$

So we have a advection component with velocity $\nu' \nabla P_f$ and a diffusion component with a coefficient $\nu' P_f$.

As part of the production $P_f$ we disintegrate older sediment at a fixed rate.

## Test 1: production transport

Suppose we have an incline in one direction, as per usual on a coastal slice. Production is happening in a circular patch in our box, with constant rate. In addition, we'll release the top 1m of sediment for further transport.

```@raw html
<details><summary>Test model</summary>
```

``` {.julia file=examples/transport/active-layer.jl}
module ActiveLayer

using Unitful
using CarboKitten.Config: Box, axes
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Utility: in_units_of
using CarboKitten.Transport.ActiveLayer: pde_stencil, Amount, Rate

<<example-active-layer>>

end
```

```@raw html
</details>
```

Our input structure facilitates a single facies, specifying an initial bedrock elevation, sediment layer and a function for a location dependent constant production rate. The transport is parametrized by a disintegration rate and a diffusion coefficient.

``` {.julia #example-active-layer}
@kwdef struct Input
    box
    Δt::typeof(1.0u"Myr")
    t_end::typeof(1.0u"Myr")
    initial_topography   # function (x::u"m", y::u"m") -> u"m"
    initial_sediment    # function (x::u"m", y::u"m") -> u"m"
    production          # function (x::u"m", y::u"m") -> u"m/s"
    disintegration_rate::typeof(1.0u"m/Myr")
    subsidence_rate::typeof(1.0u"m/Myr")
    diffusion_coefficient::typeof(1.0u"m/yr")
end
```

### Production patch
Establish a grid of 100x50, 15km on each side, dropping from 0 to 50m depth. Keeping the disintegration rate to a similar value as the production rate seems a sensible choice.

``` {.julia #example-active-layer}
production_patch(center, radius, rate) = function(x, y)
    (pcx, pcy) = center
    (x - pcx)^2 + (y - pcy)^2 < radius^2 ?
        rate :
        0.0u"m/Myr"
end

const input = Input(
    box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0u"m"),
    Δt=0.001u"Myr",
    t_end=1.0u"Myr",

    initial_topography = (x, y) -> -x / 300.0,
    initial_sediment = (x, y) -> 0.0u"m",

    production = production_patch(
        (5000.0u"m", 3750.0u"m"),
        2.0u"km",
        50.0u"m/Myr"),

    disintegration_rate = 50.0u"m/Myr",
    subsidence_rate = 50.0u"m/Myr",

    diffusion_coefficient = 10.0u"m/yr"
)
```

![Production patch on an inclining bedrock](fig/active-layer-production-patch.png)

```@raw html
<details><summary>Plotting code</summary>
```

``` {.julia .task file=examples/transport/active-layer-plot-production.jl}
#| requires: examples/transport/active-layer.jl
#| creates: docs/src/_fig/active-layer-production-patch.png
#| collect: figures

include("active-layer.jl")
using Unitful
using CarboKitten.Config: axes
using CarboKitten.Utility: in_units_of
using CairoMakie
using .ActiveLayer: input

function main()
  (x, y) = axes(input.box)
  η = input.initial_topography.(x, y')
  p = input.production.(x, y')

  fig = Figure()
  ax = Axis3(fig[1,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
  surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"), color = p |> in_units_of(u"m/Myr"))
  save("docs/src/_fig/active-layer-production-patch.png", fig)
end

main()
```

```@raw html
</details>
```

### Solving the PDE

Just as a reminder:

$$\partial_t \eta_f = \nu' \nabla P_f(x) \cdot \nabla \eta(x) + \nu' P_f(x) \nabla^2 \eta(x) + P_f(x)$$

Below is the kernel encoding a central differencing scheme i.e. `[-1, 0, 1]/(2Δx)` for first derivative and `[0 -1 0; -1 4 -1; 0 -1 0]/Δx^2` for the laplacian.

``` {.julia file=src/Transport/ActiveLayer.jl}
module ActiveLayer

using Unitful
using StaticArrays
using ...BoundaryTrait
using ...Boxes: Box
using ...Stencil: stencil!

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")

function pde_stencil(box::Box{BT}, Δt, ν, out, η, C) where {BT<:Boundary{2}}
    Δx = box.phys_scale
    d = ν * Δt

    stencil!(BT, Size(3, 3), out, η, C) do η, C
        adv = d * ((η[3, 2] - η[1, 2]) * (C[3, 2] - C[1, 2]) +
                   (η[2, 3] - η[2, 1]) * (C[2, 3] - C[2, 1])) /
              (2Δx)^2

        dif = d * C[2, 2] * (η[3, 2] + η[2, 3] + η[1, 2] +
                             η[2, 1] - 4 * η[2, 2]) / (Δx)^2

        prd = C[2, 2]

        max(0.0u"m", adv + dif + prd)
    end
end

end
```

### Model loop

Every iteration we determine the maximum disintegrated sediment. If the total amount of sediment is smaller than the maximum, then that amount is disintegrated instead. We compute the concentrations in the active layer in terms of amounts of sediment, so $P \Delta t$. Since $P$ appears in every term of the PDE, we're free to do so.

``` {.julia #example-active-layer}
mutable struct State
    time::typeof(1.0u"Myr")
    sediment::Matrix{typeof(1.0u"m")}
end

function initial_state(input)
    x, y = axes(input.box)
    State(0.0u"Myr", input.initial_sediment.(x, y'))
end

struct Frame
    t::typeof(1.0u"Myr")
    δ::Matrix{Amount}
end

function propagator(input)
    δ = Matrix{Amount}(undef, input.box.grid_size...)
    x, y = axes(input.box)
    μ0 = input.initial_topography.(x, y')
    box = input.box
    Δt = input.Δt
    disintegration_rate = input.disintegration_rate
    production = input.production
    d = input.diffusion_coefficient

    function active_layer(state)
        max_amount = disintegration_rate * Δt
        amount = min.(max_amount, state.sediment)
        state.sediment .-= amount

        production.(x, y') * Δt .+ amount
    end

    function (state)
        p = active_layer(state)
        pde_stencil(box, Δt, d, δ, state.sediment .+ μ0, p)
        return Frame(state.time, δ)
    end
end

function run_model(input)
    state = initial_state(input)
    prop = propagator(input)

    Channel{State}() do ch
        while state.time < input.t_end
            Δ = prop(state)
            state.sediment .+= Δ.δ
            state.time += input.Δt
            put!(ch, state)
        end
    end
end
```

### Running the model

We run the model with 1000 time steps but only inspect one in every 100.

![Active layer test](fig/active-layer-test.png)

```@raw html
<details><summary>Plotting code</summary>
```

``` {.julia .task file=examples/transport/active-layer-plot-result.jl}
#| requires: examples/transport/active-layer.jl
#| creates: docs/src/_fig/active-layer-test.png
#| collect: figures

module ActiveLayerPlot

include("active-layer.jl")
using CairoMakie
using Unitful
using CarboKitten.Config: axes
using CarboKitten.Utility: in_units_of
using .ActiveLayer: input, run_model

function main()
  result = Iterators.map(deepcopy,
      Iterators.filter(x -> mod(x[1], 100) == 0, enumerate(run_model(input)))) |> collect

    (x, y) = axes(input.box)
    η = input.initial_topography.(x, y') .+ result[10][2].sediment .- input.subsidence_rate * result[10][2].time
    # p = input.production.(x, y')

    fig = Figure(size=(800, 1000))
    ax = Axis3(fig[1:2,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
    surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"))

    ax2 = Axis(fig[3,1], xlabel="x (km)", ylabel="η (m)")

    for i in 1:10
        η = input.initial_topography.(x, y') .+ result[i][2].sediment .- input.subsidence_rate * result[i][2].time

        lines!(ax2, x |> in_units_of(u"km"), η[:, 25] |> in_units_of(u"m"))
    end

    save("docs/src/_fig/active-layer-test.png", fig)
end

end

ActiveLayerPlot.main()
```

```@raw html
</details>
```

Note in the bottom figure, due to sedimentation not keeping up with subsidence, the lines go down in time. We see the sediment transport being favoured to downslope areas, which is what we want. This effect could be made more extreme by increasing the disintegration rate.

## One-dimensional tests

``` {.julia file=examples/transport/plot-1d-evolution.jl}
<<plot-1d-evolution>>
```

``` {.julia #plot-1d-evolution}
using Printf: @sprintf
using Unitful: ustrip

function plot_1d_evolution!(ax::Axis, input, every=100)
	y_idx = 1
	(x, y) = box_axes(input.box)
	state = ALCAP.initial_state(input)

	plot_state() = begin
		t = state.step * input.time.Δt
		η = input.initial_topography.(x, y') .+ 
            state.sediment_height .-
            input.subsidence_rate * t
		lines!(ax, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"),
               label=@sprintf("%.3f Myr", ustrip(t)))
	end

	plot_state()
	run_model(Model{ALCAP}, input, state) do i, _
		if mod(i, every) == 0
			plot_state()
		end
	end
end

function plot_1d_evolution(input, every=100)
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

    plot_1d_evolution!(ax, input, every)

	Legend(fig[1, 2], ax)

	fig
end
```

``` {.julia file=src/Testing.jl}
module Testing

using ..CarboKitten
using ..Models.ALCAP
using Unitful
using GeometryBasics

transport_test_input(;
	initial_topography = (x, y) -> 0.0u"m",
	initial_sediment = (x, y) -> 0.0u"m",
	disintegration_rate = 50.0u"m/Myr",
	subsidence_rate = 0.0u"m/Myr",
	diffusion_coefficient = 0.0u"m/yr",
	wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr"),
    intertidal_zone = 0.0u"m") =

	ALCAP.Input(
		box = CarboKitten.Box{Coast}(grid_size=(120, 1), phys_scale=125.0u"m"),
		time = TimeProperties(
			Δt = 0.001u"Myr",
			steps = 1000),
		facies = [ALCAP.Facies(
			initial_sediment = initial_sediment,
			diffusion_coefficient = diffusion_coefficient,
			wave_velocity = wave_velocity,
			maximum_growth_rate = 0.0u"m/Myr",
			extinction_coefficient = 0.8u"m^-1",
			saturation_intensity = 60u"W/m^2"
		)],
		disintegration_rate = disintegration_rate,
		initial_topography = initial_topography,
		insolation = 400.0u"W/m^2",
		sediment_buffer_size = 5,
		depositional_resolution = 1000.0u"m",
		transport_solver = Val{:forward_euler},
        subsidence_rate = subsidence_rate,
        intertidal_zone = intertidal_zone)

end
```

### Erosion

The erosion scenario tests that sharp sediment profiles erode under diffusive transport.

![Erosion](fig/1d-erosion.svg)

``` {.julia .task file=examples/transport/1d-erosion.jl}
#| creates: docs/src/_fig/1d-erosion.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

	function initial_sediment(x, y)
	  if x < 5.0u"km"
	    return 30.0u"m"
	  end

	  if x > 10.0u"km" && x < 11.0u"km"
	    return 20.0u"m"
	  end

	  return 5.0u"m"
	end

    function main()
        CairoMakie.activate!()
        input = transport_test_input(
            initial_topography = (x, y) -> -30.0u"m",
            initial_sediment = initial_sediment,
            diffusion_coefficient = 10.0u"m/yr")

        fig = plot_1d_evolution(input, 250)
        save("docs/src/_fig/1d-erosion.svg", fig)
    end
end

Script.main()
```

### Advection

To test advective properties, we set disintegration rate to infinity, This way, surface gradients are zero and we get pure advection.

![advection](fig/1d-advection.svg)

``` {.julia .task file=examples/transport/1d-advection.jl}
#| creates: docs/src/_fig/1d-advection.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

    function gaussian_initial_sediment(x, y)
        exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
    end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    function main()
        CairoMakie.activate!()
        input = transport_test_input(
            initial_topography = (x, y)  -> -30.0u"m",
            initial_sediment = gaussian_initial_sediment,
            disintegration_rate = 50000.0u"m/Myr",
            wave_velocity = v_const(-5u"km/Myr")
        )

        fig = plot_1d_evolution(input, 250)
        save("docs/src/_fig/1d-advection.svg", fig)
    end
end

Script.main()
```

### Onshore transport

This test shows the difference between having a velocity profile without and with shear.

![onshore](fig/1d-onshore.svg)

``` {.julia .task file=examples/transport/1d-onshore.jl}
#| creates: docs/src/_fig/1d-onshore.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

    v_prof(v_max, max_depth, w) = 
        let k = sqrt(0.5) / max_depth,
            A = 3.331 * v_max,
            α = tanh(k * w),
            β = exp(-k * w)
            (A * α * β, -A * k * β * (1 - α - α^2))
        end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    v_prof_par(v_max, max_depth) = w -> let (v, s) = v_prof(v_max, max_depth, w)
        (Vec2(v, 0.0u"m/yr"), Vec2(s, 0.0u"1/yr"))
    end

    function main()
        CairoMakie.activate!()

        fig = Figure(size=(1000, 500))
        ax1 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")
        input1 = transport_test_input(
            initial_topography = (x, y) -> -x / 375.0 - 10u"m",
            initial_sediment = 10.0u"m",
            disintegration_rate = 50.0u"m/Myr",
            diffusion_coefficient = 5.0u"m/yr",
            wave_velocity = v_const(-0.5u"m/yr")
        )
        plot_1d_evolution!(ax1, input1, 250)

        ax2 = Axis(fig[1, 2], xlabel="x (km)", ylabel="η (m)")
        input2 = transport_test_input(
            initial_topography = (x, y) -> -x / 375.0 - 10u"m",
            initial_sediment = 10.0u"m",
            disintegration_rate = 50.0u"m/Myr",
            diffusion_coefficient = 5.0u"m/yr",
            wave_velocity = v_prof_par(-0.5u"m/yr", 20u"m")
        )
        plot_1d_evolution!(ax2, input2, 250)

        Legend(fig[1, 3], ax2)
        save("docs/src/_fig/1d-onshore.svg", fig)
    end
end

Script.main()
```

## Active Layer Component

```component-dag
CarboKitten.Components.ActiveLayer
```

``` {.julia file=src/Components/ActiveLayer.jl}
@compose module ActiveLayer
@mixin WaterDepth, FaciesBase, SedimentBuffer

export disintegrator, transporter, precipitation_factor

using ..Common
using CarboKitten.Transport.Advection: transport, advection_coef!, transport_dC!, max_dt
using CarboKitten.Transport.Solvers: runge_kutta_4, forward_euler
using Unitful
using GeometryBasics

@kwdef struct Facies <: AbstractFacies
    diffusion_coefficient::typeof(1.0u"m/yr") = 0.0u"m/Myr"
    wave_velocity = _ -> (Vec2(0.0u"m/Myr", 0.0u"m/Myr"), Vec2(0.0u"1/Myr", 0.0u"1/Myr"))
end

@kwdef mutable struct State <: AbstractState
    active_layer::Array{Amount, 3}
end

@kwdef struct Input <: AbstractInput
    intertidal_zone::Height = 0.0u"m"
    disintegration_rate::Rate = 50.0u"m/Myr"
    precipitation_time::Union{typeof(1.0u"Myr"), Nothing} = nothing
    transport_solver = Val{:RK4}
    transport_substeps = :adaptive 
end

courant_max(::Type{Val{:RK4}}) = 2.0
courant_max(::Type{Val{:forward_euler}}) = 1.0

transport_solver(f, _) = f
transport_solver(::Type{Val{:RK4}}, box) = runge_kutta_4(typeof(1.0u"m"), box)
transport_solver(::Type{Val{:forward_euler}}, _) = forward_euler

function precipitation_factor(input::AbstractInput)
    if input.precipitation_time === nothing
        return 1.0
    else
        return 1.0 - exp(input.time.Δt * log(1/2) / input.precipitation_time)
    end
end

function adaptive_transporter(input)
    solver = transport_solver(input.transport_solver, input.box)

    w = water_depth(input)
    box = input.box
    Δt = input.time.Δt  # / input.transport_substeps
    fs = input.facies
    adv = Matrix{Vec2{Rate}}(undef, box.grid_size...)
    rct = Matrix{typeof(1.0u"1/Myr")}(undef, box.grid_size...)
    dC = Matrix{Rate}(undef, box.grid_size...)
    cm = courant_max(input.transport_solver)
    iz = input.intertidal_zone

    return function (state)
        wd = w(state)
        wd .+= iz

        C = state.active_layer
        for (i, f) in pairs(fs)
            advection_coef!(box, f.diffusion_coefficient, f.wave_velocity, wd, adv, rct)
            m = max_dt(adv, box.phys_scale, cm)
            steps = ceil(Int, Δt / m)

            # @debug "step $(state.step) - transport substeps $(steps)"
            subdt = Δt / steps
            for j in 1:steps
                solver(
                    (C, _) -> transport_dC!(input.box, adv, rct, C, dC),
                    view(C, i, :, :), TimeIntegration.time(input, state), subdt)
            end
        end

        for i in eachindex(C)
            if C[i] < zero(Amount)
                C[i] = zero(Amount)
            end
        end
    end
end

"""
    disintegrator(input) -> f!

Prepares the disintegration step. Returns a function `f!(state::State)`. The returned function
modifies the state, popping sediment from the `sediment_buffer` and returns an array of `Amount`.
"""
function disintegrator(input)
    max_h = input.disintegration_rate * input.time.Δt
    w = water_depth(input)
    output = Array{Float64,3}(undef, n_facies(input), input.box.grid_size...)
    depositional_resolution = input.depositional_resolution
    iz = input.intertidal_zone

    return function (state)
        wn = w(state)
        wn .+= iz
        h = min.(max_h, state.sediment_height)
        h[wn.<=0.0u"m"] .= 0.0u"m"

        @assert all(h .<= max_h)
        state.sediment_height .-= h
        pop_sediment!(state.sediment_buffer, h ./ depositional_resolution .|> NoUnits, output)
        return output .* depositional_resolution
    end
end

"""
    transporter(input::Input) -> f

Prepares the transportation step. Returns a function `f(state::State, active_layer)`,
transporting the active layer, returning a transported `Amount` of sediment.
"""
function transporter(input)
    if input.transport_substeps == :adaptive
        return adaptive_transporter(input)
    end

    solver = transport_solver(input.transport_solver, input.box)

    w = water_depth(input)
    box = input.box
    Δt = input.time.Δt / input.transport_substeps
    steps = input.transport_substeps
    fs = input.facies
    iz = input.intertidal_zone

    return function (state)
        wd = w(state)
        wd .+= iz

        C = state.active_layer
        for (i, f) in pairs(fs)
            for j in 1:steps
                solver(
                    (C, _) -> transport(
                        input.box, f.diffusion_coefficient, f.wave_velocity,
C, wd),
                    view(C, i, :, :), TimeIntegration.time(input, state), Δt)
            end
        end

        for i in eachindex(C)
            if C[i] < zero(Amount)
                C[i] = zero(Amount)
            end
        end
    end
end

end
```

## Intertidal zone

You may pass the `intertidal_zone` argument to define a height above sea-level where sediment is still transported.

The following test has three sediment peaks, one below sea-level, one in the intertidal zone and one above the intertidal zone. In the latter two cases, sediment is transported.

The `transport_test_input` has a box of $120 \times 1$ with a resolution of 125m (15km total).
We divide the domain into three sections: dry ($h > 10m$), intertidal ($10m > h > 0m$) and wet zones ($0m > h$). The dry zone never has transport, the wet zone always has transport, but the intertidal zone only has transport enabled when we set `intertidal_zone` to `10u"m"`.

![intertidal zone test](fig/1d-intertidal.svg)

``` {.julia .task file=examples/transport/1d-intertidal-zone.jl}
#| creates: docs/src/_fig/1d-intertidal.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

    function three_peaks(x, y)
        sum(exp(-(x-μ)^2 / (2 * (0.5u"km")^2)) * 9.0u"m"
            for μ in [2.5u"km", 7.5u"km", 12.5u"km"])
    end

    function staircase(dx, dy, y0)
        return function (x, _)
            floor(x / dx) * dy + y0
        end
    end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    function main()
        CairoMakie.activate!()

        input = transport_test_input(
            initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
            initial_sediment = three_peaks,
            disintegration_rate = 50000.0u"m/Myr",
            wave_velocity = v_const(-1u"km/Myr"),
            intertidal_zone = 10u"m"
        )

        fig = plot_1d_evolution(input, 500)
        save("docs/src/_fig/1d-intertidal.svg", fig)
    end
end

Script.main()
```

The following tests that we see the expected behaviours both without an intertidal zone (`input1`) and with (`input1`). The regions split at grid locations 40 and 80, so we test slices `10:30`, `50:70` and `90:110`. This ommits the boundaries where numeric artifacts could be encountered.

``` {.julia file=test/Transport/IntertidalZoneSpec.jl}
@testset "CarboKitten.Transport.IntertidalZone" begin
    using CarboKitten
    using CarboKitten.Testing: transport_test_input

    function end_sediment_height(input)
        state = ALCAP.initial_state(input)
        run_model((_, _) -> (), Model{ALCAP}, input, state)
        return state.sediment_height
    end

    function three_peaks(x, y)
        sum(exp(-(x-μ)^2 / (2 * (0.5u"km")^2)) * 9.0u"m"
            for μ in [2.5u"km", 7.5u"km", 12.5u"km"])
    end

    function staircase(dx, dy, y0)
        return function (x, _)
            floor(x / dx) * dy + y0
        end
    end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    input1 = transport_test_input(
        initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
        initial_sediment = three_peaks,
        disintegration_rate = 50000.0u"m/Myr",
        wave_velocity = v_const(-1u"km/Myr"),
        intertidal_zone = 0u"m"
    )

    output1 = end_sediment_height(input1)[:, 1]

    input2 = transport_test_input(
        initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
        initial_sediment = three_peaks,
        disintegration_rate = 50000.0u"m/Myr",
        wave_velocity = v_const(-1u"km/Myr"),
        intertidal_zone = 10u"m"
    )

    output2 = end_sediment_height(input2)[:, 1]

    @test output1[10:30] ≈ output1[50:70] atol=0.01u"m"
    @test !isapprox(output1[50:70], output1[90:110], atol=1.0u"m")

    @test !isapprox(output2[10:30], output2[50:70], atol=1.0u"m")
    @test output2[50:70] ≈ output2[90:110] atol=0.01u"m"
end
```
