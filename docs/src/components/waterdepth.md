# Water Depth

```component-dag
CarboKitten.Components.WaterDepth
```

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

## Input

- `initial_topography(x, y)` (a.k.a. initial depth) should be a function taking two coordinates in units of meters, returning an elevation also in meters.
- `sea_level(t)` should be a function taking a time in millions of years (Myr) returning the eustatic sealevel. This could also be an interpolated table.
- `subsidence_rate` the rate of subsidence. This is one of:
    - a constant rate, e.g. `50.0u"m/Myr"`;
    - a `TimeVariableSubsidenceRate(t -> ...)`, wrapping a function of time returning a rate;
    - a `SpaceTimeVariableSubsidenceRate((x, y, t) -> ...)`, wrapping a function of position and time returning a rate.

The signs of these quantities should be such that the following equation holds:

$$T + E = S + W,$$

saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

### Variable subsidence

The subsidence rate is integrated with a forward Euler step: at step $k$ the bathymetry is lowered by $\sigma(x, y, t_k) \Delta t$, where $t_k = t_0 + k \Delta t$. For instance, to model a margin that tilts seaward, where fault activity dies out over the course of the run, we could write:

```julia
subsidence_rate = SpaceTimeVariableSubsidenceRate(
    (x, y, t) -> 20.0u"m/Myr" + x * 5.0u"1/Myr" * exp(-t / 0.5u"Myr"))
```

Note that the function is evaluated at the cell centers, and that, like with `initial_topography`, the function needs to be vectorizable using the dot syntax.

Since the subsidence is no longer a single number, post-processing can't simply multiply the rate with the duration of the run. Instead, the `cumulative_subsidence` function computes the total amount of subsidence for each grid cell, by running the `subsider` on an empty state. This guarantees that the result agrees exactly with what the model does. The output header stores this total, so that `bedrock(header)` gives the correct bedrock for any type of subsidence.

``` {.julia file=src/Components/WaterDepth.jl}
@component module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
import ..Common: write_header
using HDF5
using ..TimeIntegration: time, time_axis

export water_depth, subsider, initial_topography, cumulative_subsidence
export TimeVariableSubsidenceRate, SpaceTimeVariableSubsidenceRate

@kwdef struct Input <: InputPart
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    subsidence_rate = 0.0u"m/Myr"
end

@kwdef mutable struct State <: StatePart
    bathymetry::Matrix{Height}
end

function initial_state(input::AbstractInput)
    bathymetry = initial_topography(input)
    return record_type(:State, WaterDepth)(step=0, bathymetry=bathymetry)
end

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

struct TimeVariableSubsidenceRate
    f
end

struct SpaceTimeVariableSubsidenceRate
    f
end

subsider(input::AbstractInput) = subsider(input, input.subsidence_rate)

function subsider(input::AbstractInput, s::SpaceTimeVariableSubsidenceRate)
    get_time = time(input)
    x, y = box_axes(input.box)
    dt = input.time.Δt

    return function (state::AbstractState)
        t = get_time(state)
        Δσ = s.f.(x, y', t) * dt
        state.bathymetry .-= Δσ
    end
end

function subsider(input::AbstractInput, s::TimeVariableSubsidenceRate)
    get_time = time(input)
    dt = input.time.Δt

    return function (state::AbstractState)
        t = get_time(state)
        Δσ = s.f(t) * dt
        state.bathymetry .-= Δσ
    end
end

function subsider(input::AbstractInput, s::Quantity)
    Δσ = input.subsidence_rate * input.time.Δt

    return function (state::AbstractState)
        state.bathymetry .-= Δσ
    end
end

"""
    cumulative_subsidence(input)

Total amount of subsidence over the entire run, for each grid cell.
"""
function cumulative_subsidence(input::AbstractInput)
    state = record_type(:State, WaterDepth)(step=0, bathymetry=zeros(Amount, input.box.grid_size...))
    subside! = subsider(input)
    for _ in 1:input.time.steps
        subside!(state)
        state.step += 1
    end
    return .-state.bathymetry
end

function water_depth(input::AbstractInput)
    sea_level = input.sea_level
    get_time = time(input)

    return function (state::AbstractState)
        t = get_time(state)
        return sea_level(t) .- state.bathymetry
    end
end

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))
    set_attribute(output, "total_subsidence", cumulative_subsidence(input) |> in_units_of(u"m"))
    if input.subsidence_rate isa Quantity
        set_attribute(output, "subsidence_rate", input.subsidence_rate |> in_units_of(u"m/Myr"))
    end
end

end
```

## Tests

We test that each type of subsidence rate gives the expected bathymetry after a few steps, and that `cumulative_subsidence` agrees with stepping the `subsider`.

``` {.julia file=test/Components/WaterDepthSpec.jl}
using CarboKitten
import CarboKitten.Components.WaterDepth as WD
using CarboKitten.Components.Common: record_type
using CarboKitten.Components.WaterDepth:
    TimeVariableSubsidenceRate, SpaceTimeVariableSubsidenceRate

@testset "Components/WaterDepth" begin
    input = record_type(:Input, WD)(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=10),
        sea_level = t -> 2.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 5.0u"m/Myr"
    )
    state = WD.initial_state(input)


    @test all(state.bathymetry .== WD.initial_topography(input))

    sub! = WD.subsider(input)
    sub!(state)
    @test all(state.bathymetry.==-15.0u"m")

    wd = WD.water_depth(input)
    @test all(wd(state) .== 17.0u"m")

    @test all(WD.cumulative_subsidence(input) .≈ 50.0u"m")
end

@testset "Components/WaterDepth/variable-subsidence" begin
    box = Box{Periodic{2}}(grid_size=(10, 5), phys_scale=1.0u"m")
    time = TimeProperties(t0=1.0u"Myr", Δt=1.0u"Myr", steps=4)
    x, y = box_axes(box)

    function run_subsider(input)
        state = record_type(:State, WD)(step=0, bathymetry=zeros(typeof(1.0u"m"), box.grid_size...))
        sub! = WD.subsider(input)
        for _ in 1:input.time.steps
            sub!(state)
            state.step += 1
        end
        return state.bathymetry
    end

    @testset "time variable" begin
        # rate is 1, 2, 3, 4 m/Myr at t = 1, 2, 3, 4 Myr
        input = record_type(:Input, WD)(
            box = box, time = time,
            subsidence_rate = TimeVariableSubsidenceRate(t -> t * 1.0u"m/Myr^2"))
        @test all(run_subsider(input) .≈ -10.0u"m")
        @test all(WD.cumulative_subsidence(input) .≈ 10.0u"m")
    end

    @testset "space-time variable" begin
        f(x, y, t) = x * t * 1.0u"1/Myr^2" + y * 1.0u"1/Myr"
        input = record_type(:Input, WD)(
            box = box, time = time,
            subsidence_rate = SpaceTimeVariableSubsidenceRate(f))
        expected = x .* 10.0u"Myr^-1" .* 1.0u"Myr" .+ y' .* 4.0u"Myr^-1" .* 1.0u"Myr"
        @test run_subsider(input) ≈ .-expected
        @test WD.cumulative_subsidence(input) ≈ expected
    end

    @testset "constant rate is a special case" begin
        constant = record_type(:Input, WD)(box = box, time = time, subsidence_rate = 3.0u"m/Myr")
        timevar = record_type(:Input, WD)(box = box, time = time,
            subsidence_rate = TimeVariableSubsidenceRate(_ -> 3.0u"m/Myr"))
        spacetimevar = record_type(:Input, WD)(box = box, time = time,
            subsidence_rate = SpaceTimeVariableSubsidenceRate((_, _, _) -> 3.0u"m/Myr"))
        @test run_subsider(timevar) ≈ run_subsider(constant)
        @test run_subsider(spacetimevar) ≈ run_subsider(constant)
    end
end
```
