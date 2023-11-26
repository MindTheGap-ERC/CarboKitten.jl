# Submarine Transport

## Input parameters

### Facies properties
Facies now have an additional property: average grain size.

``` {.julia #cat-facies}
struct Facies
  viability_range::Tuple{Int,Int}
  activation_range::Tuple{Int,Int}

  maximum_growth_rate::Float64
  extinction_coefficient::Float64
  saturation_intensity::Float64

  grain_size::Float64
end
```

### Buffer properties
The sediment output will be stored in a buffer with fixed spatial resolution, typically a few pixels per meter depth. This buffer needs to store percentages for each facies, and age, so size will be $x\ \times\ y\ \times\ z\ \times\ (n_{\rm spec} + 1) \ \times {\rm f32}$, typically a few tens of MB. New properties are `Δz` for the buffer resolution and `buffer_depth` for the size of the buffer.

``` {.julia #cat-input}
@kwdef struct Input
  sea_level
  subsidence_rate
  initial_depth

  grid_size::NTuple{2,Int}
  phys_scale::Float64
  Δt::Float64

  time_steps::Int
  write_interval::Int

  facies::Vector{Facies}
  insolation::Float64

  Δz::Float64
  buffer_depth::Int
end
```

## Model structure
The addition of sediment transport to the model requires an update to the simple structure we had when only production was involved. Remember, we have a **state** $S$, a **propagator**,

$$P_i: S \to \Delta,$$

and a function $U$ that *updates* the state with the given frame,

$$U: (S, \Delta) \to S.$$

In the new setup, our state will also contain a layer of history upon which we can erode and transport sediment. Also, each time step will consist of two deltas, one for production and one for transport. (If you want to disperse production before sedimentation, this is still considered production stage).

``` {.julia #cat-frame}
struct ProductFrame
  production::Array{Float64,3}        # x y f
end
```

``` {.julia #cat-frame}
struct TransportFrame
  disintegration::Array{Float64,2}   # x y
  deposition::Array{Float64,3}       # x y f
end
```

### Buffer
There are several choices on how to structure the sediment buffer. We can grow sediment from the bottom of the buffer, but that requires keeping an array of pointers to the bottom depth where new material is deposited.

Another choice is to keep the sea floor in the same layer of the buffer, and copy down sediment when a layer is full. We may need to implement both to see which is more efficient. 

``` {.julia file=test/SedimentStackSpec.jl}
@testset "SedimentStack" begin
  stack = zeros(Float64, 10, 3)
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
end
```

``` {.julia file=src/SedimentStack.jl}
module SedimentStack

export push_sediment!, pop_sediment!

function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
  Δ = sum(parcel)
  bucket = sum(col[1, :])

  if bucket + Δ < 1.0
    col[1,:] .+= parcel
    return
  end

  frac = parcel ./ Δ
  col[1,:] .+= frac .* (1.0 - bucket)
  Δ -= (1.0 - bucket)
  n = floor(Int64, Δ)

  col[n+2:end,:] = col[1:end-n-1,:]
  na = [CartesianIndex()]
  col[2:n+1,:] .= frac[na,:]
  Δ -= n

  col[1,:] .= frac .* Δ
end

@inline function pop_fraction(col::AbstractMatrix{F}, Δ::F) where F <: Real
  bucket = sum(col[1,:])
  @assert Δ < bucket "pop_fraction can only pop from the top cell"
  parcel = (Δ / bucket) .* col[1,:]
  col[1,:] .-= parcel
  return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
  bucket = sum(col[1,:])
  if Δ < bucket
    return pop_fraction(col, Δ)
  end

  parcel = copy(col[1,:])
  Δ -= bucket
  n = floor(Int64, Δ)

  parcel .+= sum(col[2:n+1,:]; dims=1)'
  col[1:end-n-1, :] = col[n+2:end, :]
  col[end-n-1:end, :] .= 0
  Δ -= n

  parcel .+= pop_fraction(col, Δ)
  return parcel
end

end # module
```

#### Current choice
Implementation will be such that depth $z$ is the smallest axis, such that depth transects are contiguous in memory. The sediment buffer will be four dimensional containing fractions for each species, i.e. $\in [0, 1]$, and one additional slice to record the time.

``` {.julia #cat-state}
mutable struct State
  time::Float64
  height::Array{Float64,2}            # x y
  sediment::Array{Float64,4}          # x y z (f..., t)
end
```

### Delta
Each time step we compute the production as before, resulting in a $\Delta_{\prod}$. What has changed, is how we update the state from this propagator.

``` {.julia #cat-update}
function deposit_material(
    grid_size::NTuple{2, Int},
    Δt::Float64,
    Δz::Float64,
    s::State,
    facies::T) where T <: AbstractArray{Float64,3}

  Threads.@threads for idx in CartesianIndices(grid_size)
    prod = facies[Tuple(idx)..., :]
    Δh = sum(prod) .* Δt
    fractions = prod ./ sum(prod)
    column = s.sediment[Tuple(idx)..., :, :]
    bucket = sum(column[1, :])
    if bucket + Δh / Δz > 1.0
      column[1, :] .+= fractions .* (1.0 - bucket)
      column[2:end, :] = column[1:end-1, :]
      column[1, :] = fractions .* (Δh / Δz - 1.0)
    else
      column[1, :] .+= fractions .* (Δh / Δz)
    end

    s.height[idx] .-= Δh
  end
end

function production_updater(input::Input)
  function (s::State, Δ::ProductFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.production)
  end
end
```

``` {.julia #cat-update}
function transport_updater(input::Input)
  function (s::State, Δ::TransportFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.deposition)
  end
end
```

``` {.julia #cat-update}
function time_updater(input::Input)
  function (s::State, Δ::Nothing = nothing)
    s.height .+= input.subsidence_rate * input.Δt
    s.time += input.Δt
  end
end
```

``` {.julia file=src/CATP.jl}
module CATP

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.Stencil
using CarboKitten.Burgess2013.CA

using HDF5

using .Iterators: drop, peel, partition, map, take
<<cat-facies>>
<<cat-input>>
<<cat-state>>
<<cat-frame>>
<<cat-update>>
<<ca-prod-propagator>>

end  # CATP
```

## Science description
In this model, sediment is entrained (lifted from the surface and transported) if the shear stress $\tau$ exceeds the critical shear stress $\tau_{e}$.

The shear stress acting on the sediment $\tau_{total}$ can be simplified to a combination of shear stress from currents, waves and slopes:

$\tau_{total} = \tau_{current} + \tau_{wave} + \tau_{slope}$

### Shear stress due to slope

$\tau_{slope} = \Delta \rho g D \sin( \alpha ) ( \underline{G} )/G$

$\Delta \rho$ is excess density of the sediment, *g* is gravitational acceleration, *D* is grain size and $\alpha$ the max slope angle parallel to the gradient vector $\underline{G}$.

![Eq 17](https://github.com/MindTheGap-ERC/CarboKitten/assets/18270548/47d9ba18-0431-4286-8f94-00210b987e1e)


### Shear stress due to wave and currents

This can be simplified to a unidirectional shear stress vector field where the shear stress decreases with depth in a similar way as the production rate $D_{0}$.

$\tau_{wave/current} = \tau_{wave/current}(z)$

![Fig 1](https://github.com/MindTheGap-ERC/CarboKitten/assets/18270548/30790c25-5174-4dd4-ae0e-b555eaa7c4a1)

The depth-dependence function is defined piecewise between water surface and two parameter depths $z_{t}$ and $z_{b}$.
Note that this production curve is different from the one we use from Bosscher and Schlager (1992).  

### Shear stress is different for different sediment types

Different grain sizes have different critical shear stress $\tau_{e}$. Each lithofacies corresponds to different sediment producers and that translates into different grain size distributions. For example, corals build massive skeletons that form wave-resistant reefs and they disintegrate into coarse rubble. Algae that live in lagoons disintegrate into mud.

CARBONATE 3D simplifies the range of grain sizes into two: *grains* and *matrix*, which means: medium grains (2 mm and larger) and mud (a liberal definition: grains < 2 mm). Each lithofacies will disintegrate into its specific distribution of these two fractions. For example, the reef facies would produce 70% of grains and 30% of matrix and the lagoon facies 100% matrix. This requires each grid cell to have an attribute describing its sediment composition.

For a given sediment type, critical shear stress can be approximated with:

$\tau_{e} = \Delta \rho g D \sin( \alpha_{c} )$

Where $\alpha_{c}$ is the angle of repose. There are other methods based on empirical properties of different sediment types, but they are subject to various constraints discussed in Warrlich, Waltham, and Bosence (2002) p. 384-385.

## Transport

When described as function of shear stress, transport starts when $|\tau_{total}| > \tau_{e}$ for a given sediment type and stops when $|\tau_{total}| = \tau_{e}$.

Deposition rate is proportional to sediment load *L*:

$L = L_{0}e^{-(S-S_{d}/X)}$

$L_{0}$ is the initially entrained sediment load, *X* the characteristic transport distance (user input) and $S_{d}$ is the point along the path *S* where sedimentation starts.

Further details on calculating transport direction are provided in the thesis Warrlich (2000) and we know that Georg used a separate array for calculating it, this is not documented in detail, but it's possible to ask him.

## Sediment available for transport

After deposition, sediment (loose) turns into rock (solid) through the process of diagenesis. So it is only available for transport for a limited time: either directly after it had been deposited or after it had been made available from rock through disintegration. 

CARBONATE 3D uses the following order of steps (ignore anything referring to "clastic sediments"):

![flowchart](https://github.com/MindTheGap-ERC/CarboKitten/assets/18270548/e4d2a81c-7c69-432f-8371-fecd308c934a)


It is not clear to me why entrainment of disintegrated sediment precedes the production of new sediment. One might consider that disintegration acts upon rock exposed after loose sediment is entrained completely. So if in any grid cell sediment is produced and then completely removed, then this grid cell can be subject to disintegration and further sediment removal. Something like this:

1. Production leads to new loose sediment
2. Loose sediment is transported
2. In grid cells where no production took place or in which all new sediment has been removed, disintegration occurs, which also leads to new loose sediment
3. Loose sediment is transported

That also means transport calculated twice. Perhaps we could ask Georg about it.

### Disintegration

Only disintegration through hydrodynamic (wave and current) and biological activity is included. The thesis by Warrlich (2000) included disintegration as a result of oversteepening of depositional slopes, but this has been abandoned in the Warrlich, Waltham, and Bosence (2002) paper. Possibly because sediment transport in each step prevents from oversteepened sloped forming in the first place. We follow the paper and ignore it here.

Disintegration rate *B* is modeled similarly to production rate:

![disintegration](https://github.com/MindTheGap-ERC/CarboKitten/assets/18270548/0d50309b-fea9-48c0-8831-db8fc1835f87)
The article by Warrlich, Waltham, and Bosence (2002) includes a correcting factor $U_{E}(x,y)$ that reflects the proximity to open water. This can probably be omitted initially.

Values of parameters listed here are in Fig. 4.4 p. 39-40 of Warrlich (2000).

This is based on following references:

Bosscher, Hemmo, and Wolfgang Schlager. 1992. “Computer Simulation of Reef Growth.” Sedimentology 39 (3): 503–12. https://doi.org/https://doi.org/10.1111/j.1365-3091.1992.tb02130.x.

Warrlich, G. M. D. 2000. “3D Computer Forward Modelling of Carbonate Platform Evolution.” {PhD}, London: Royal Holloway University of London.

Warrlich, G. M. D., D. A. Waltham, and D. W. J. Bosence. 2002. “Quantifying the Sequence Stratigraphy and Drowning Mechanisms of Atolls Using a New 3‐D Forward Stratigraphic Modelling Program (CARBONATE 3D).” Basin Research 14 (3): 379–400. https://doi.org/10.1046/j.1365-2117.2002.00181.x.

Sadly Markdown support here does not extend to BibTex ;-) And two equations are not rendered properly for some reason.


