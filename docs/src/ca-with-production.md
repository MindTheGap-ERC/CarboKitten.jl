# Combining CA with production
This model combines BS92 production with the B13 cellular automaton.

We need as an input:

- initial depth (function of space)
- sea-level curve (function of time)
- subsidence (function of time)

These should all behave as a functions, but could also be some interpolated data. The signs of these quantities should be such that the following equation holds:

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia #ca-prod-input}
struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

struct Input
    sea_level
    subsidence_rate
    initial_depth

    grid_size::NTuple{2,Int}
    phys_scale::Float64
    Δt::Float64
    time_steps::Int

    facies::Vector{Facies}
    insolation::Float64
end

function production(insolation::Float64, facies::Facies, water_depth::Float64)
    gₘ = facies.maximum_growth_rate
    I₀ = insolation
    Iₖ = facies.saturation_intensity
    w = water_depth
    k = facies.extinction_coefficient
    return gₘ * tanh(I₀/Iₖ * exp(-w * k))
end
```

Each iteration of the model, we produce a `Frame`.

``` {.julia #ca-prod-frame}
struct Frame
    production :: Array{Float64, 3}
end
```

The frame is used to update a *state* $S$. The frame should be considered a delta for the state. So, we can reproduce the height at each step from the frames.

``` {.julia #ca-prod-state}
struct State
    time :: Float64
    height :: Array{Float64, 2}
end
```

From a dynamical modeling point of view, CarboCAT operates analogous to a forward Euler integration scheme, where some components are actually a discrete model. This means we have one function that generates a `Frame` from a `State`, called the *propagator* $P$ (this is our own nomenclature),

$$P_i: S \to \Delta.$$

The suffix $i$ here is used to indicate that the propagator depends on the input.

``` {.julia #ca-prod-model}
function initial_state(input::Input)  # -> State
    <<ca-prod-init>>
end

function propagator(input::Input)
    <<ca-prod-init-propagator>>
    function (s::State)  # -> Frame
        <<ca-prod-propagate>>
    end
end
```

We'll have a second function $U$ that *updates* the state with th given frame,

$$U: (S, \Delta) \to S.$$

In practice however, the update function changes the state in-place.

``` {.julia #ca-prod-model}
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, d::Frame)
        <<ca-prod-update>>
    end
end
```

## Init

We fill the height map with the initial depth function. It is assumed that the height only depends on the second index.

``` {.julia #ca-prod-init}
# n_facies = length(input.facies)
height = zeros(Float64, input.grid_size...)
for i in CartesianIndices(height)
    height[i] = input.initial_depth(i[2] * input.phys_scale)
end
return State(height, 0.0)
```

## Propagator

``` {.julia #ca-prod-init-propagator}
n_facies = length(input.facies)
init = rand(0:n_facies, input.grid_size...)
ca = run_ca(Reflected{2}, input.facies, init, 3)

function water_depth(s::State)
    input.sea_level(s.time) - s.height
end
```

``` {.julia #ca-prod-propagate}
result = zeros(Float64, input.grid_size..., n_facies)
facies_map = take!(ca)
w = water_depth(s)
for idx in CartesianIndices(facies_map)
    if facies_map[idx] == 0
        continue
    end
    result[Tuple(idx)..., f] = production(input.insolation, input.facies[facies_map[idx]], w[idx])
end
return result
```

## Updater

``` {.julia #ca-prod-update}
for f in 1:n_facies
    s.height .+= d[f]
end
s.height .-= input.subsidence_rate * input.Δt
s.time += input.Δt
```

## Loop

``` {.julia file=examples/ca-with-prod.jl}
using CarboKitten
using CarboKitten.Burgess2013.CA
using CarboKitten.Stencil: Reflected
using CarboKitten.Utility
using CarboKitten.BS92: g

using HDF5

<<ca-prod-input>>
<<ca-prod-frame>>
<<ca-prod-state>>


cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))

function rules(facies::Vector{Facies})
    function (neighbourhood::Matrix{Int}, order::Vector{Int})
        cell_facies = neighbourhood[3, 3]
        neighbour_count(f) = sum(neighbourhood .== f)
        if cell_facies == 0
            for f in order
                n = neighbour_count(species)
                (a, b) = facies[f].activation_range
                if a <= n && n <= b
                    return species
                end
            end
            0
        else
            n = neighbour_count(cell_facies) - 1
            (a, b) = facies[cell_facies].viability_range
            (a <= n && n <= b ? cell_facies : 0)
        end
    end    
end

function run_ca(::Type{B}, facies::Vector{Facies}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    r = rules(facies)
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), r)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

<<ca-prod-model>>

function main(input::Input)
    s = initial_state(input)
    p = propagator(input)
    u = updater(input)

    for i in 1:input.time_steps
        Δ = p(s)
        u(s, Δ)
    end
end
```
