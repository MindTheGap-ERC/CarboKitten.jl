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
  excess_density::Float64
  # critical_angle::Float64

  # magnitude of critical stress, should be ΔρDg sin α where α is the critical angle
  critical_stress::Float64
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
  boundary::Boundary{2}
  phys_scale::Float64     # km / pixel
  Δt::Float64             # Ma / timestep ?

  time_steps::Int
  write_interval::Int

  facies::Vector{Facies}
  insolation::Float64

  Δz::Float64
  buffer_depth::Int
  # To estimate the average grain size and overdensity of the sediment (and
  # thereby the shear stress), we need to peek at the composition of the top
  # layer, the `top_layer_thickness` parameter controls how deep we sample.
  top_layer_thickness::Float64
  # wave shear stress function, is a function of depth, and should return a Vec2
  wave_shear_stress
  # Gravitational accelleration
  g::Float64
  # for the transport step, parcels of sediment are converted into particles
  # by setting this parameter to something >1 you can control how many particles
  # are created on each axis. So number of particles is grid width * height * transport_subsample^2.
  transport_subsample::Int
end
```

## Model structure
The addition of sediment transport to the model requires an update to the simple structure we had when only production was involved. Remember, we have a **state** $S$, a **propagator**,

$$P_i: S \to \Delta,$$

and a function $U$ that *updates* the state with the given frame,

$$U: (S, \Delta) \to S.$$

In the new setup, our state will also contain a layer of history upon which we can erode and transport sediment. Also, each time step will consist of two deltas, one for production and one for transport. (If you want to disperse production before sedimentation, this is still considered production stage).

In fact, we should have one propagator computing the magnitude of transport, an updater removing that sediment, returning also the composition of the transported sediment and then another updater depositing that sediment again.

$$P_{\rm transport}: S \to \Delta_1$$
$$U_{\rm remove}: (S, \Delta_1) \to (S, \Delta_2),\ {\rm and}\ U_{\rm resettle}: (S, \Delta_2) \to S$$

Then, we can say $$U_{\rm transport} = U_{\rm resettle} \circ U_{\rm remove}$$ and retain our nice structure.

``` {.julia #cat-frame}
struct ProductFrame
  production::Array{Float64,3}        # facies x y
end

# struct Properties
#   excess_density::Float64
#   grain_size::Float64
# end

const Particle = Transport.Particle{Nothing}
```

``` {.julia #cat-frame}
struct ReductionFrame
  disintegration::Array{Float64,2}   # x y
end
```

### Buffer
There are several choices on how to structure the sediment buffer. We can grow sediment from the bottom of the buffer, but that requires keeping an array of pointers to the bottom depth where new material is deposited.

Another choice is to keep the sea floor in the same layer of the buffer, and copy down sediment when a layer is full. We may need to implement both to see which is more efficient. 

We define two functions `push_sediment!` and `pop_sediment!`. Given a $s \times n$ matrix, where $n$ is the number of facies types and $s$ is the depth of the stack, we can grow and shrink sediment. These functions are unit-free, setting $\Delta z$ to be equal to 1.

``` {.julia file=test/SedimentStackSpec.jl}
@testset "SedimentStack" begin
  stack = zeros(Float64, 10, 3)
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
  push_sediment!(stack, [0.0, 2.0, 0.0])   # (0 0.5) (0 1) (0.5 0.5) (1 0) ...
  @test pop_sediment!(stack, 2.0) == [0.25, 1.75, 0.0]
  @test pop_sediment!(stack, 1.5) == [1.25, 0.25, 0.0]
end

@testset "SedimentArray" begin
  sediment = zeros(Float64, 10, 3, 5, 5)
  for x in 1:10
    production = rand(3, 5, 5)
    push_sediment!(sediment, production)
  end
  a = peek_sediment(sediment, 1.0)
  @test all(sum(a; dims=1) .≈ 1.0)
end
```

```@raw html
<details><summary>SedimentStack impl</summary>
```

``` {.julia file=src/SedimentStack.jl}
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment

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

function push_sediment!(sediment::AbstractArray{F, 4}, p::AbstractArray{F, 3}) where F <: Real
  _, x, y = size(p)
  @views for i in CartesianIndices((x, y))
    push_sediment!(sediment[:, :, i[1], i[2]], p[:, i[1], i[2]])
  end
end

@inline function pop_fraction(col::AbstractMatrix{F}, Δ::F) where F <: Real
  bucket = sum(col[1,:])
  @assert Δ < bucket "pop_fraction can only pop from the top cell"
  parcel = (Δ / bucket) .* col[1,:]
  col[1,:] .-= parcel
  return parcel
end

function peek_sediment(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
  bucket = sum(col[1,:])
  if Δ < bucket
    parcel = (Δ / bucket) .* col[1,:]
    return parcel 
  end

  parcel = copy(col[1,:])
  Δ -= bucket
  n = floor(Int64, Δ)

  parcel .+= sum(col[2:n+1,:]; dims=1)'
  Δ -= n

  last_bit = (Δ / sum(col[n+2,:])) .* col[n+2,:]
  parcel .+= last_bit

  return parcel
end

function peek_sediment(sediment::AbstractArray{F,4}, Δ::F) where F <: Real
  _, f, x, y = size(sediment)
  out = Array{F, 3}(undef, f, x, y)
  for i in CartesianIndices((x, y))
    out[:, i[1], i[2]] = peek_sediment(@view(sediment[:, :, i[1], i[2]]), Δ)
  end
  return out
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

```@raw html
</details>
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
    Δz::Float64,
    s::State,
    material::AbstractArray{Float64,3})

  Threads.@threads for idx in CartesianIndices(grid_size)
    production = material[Tuple(idx)..., :] .* (1.0 / Δz)
    push_sediment!(s.sediment, production)
    s.height[idx] .-= sum(material[Tuple(idx)..., :]) * Δt
  end
end

function remove_material(
    grid_size::NTuple{2, Int},
    n_facies::Int,
    Δz::Float64,
    s::State,
    thickness::AbstractArray{Float64,2})

  material = Array{Float64,3}(undef, grid_size..., n_facies)
  Threads.@threads for idx in CartesianIndices(grid_size)
    material[Tuple(idx)..., :] = pop_sediment!(s.sediment, thickness[idx])
  end
  return material
end

function deposit_updater(input::Input)
  function (s::State, Δ::ProductFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.production)
  end
end
```

``` {.julia #cat-update}
function remove_updater(input::Input)
  function (s::State, Δ::ReductionFrame)
    loose = remove_material(input.grid_size, length(input.facies), input.Δz, s, Δ.disintegration)
    # need physics to disperse loose material
    return ProductFrame(loose)
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

## Distintegration
We have two propagators: one to compute production, that goes unchanged, the other to compute **the amount of disintegration**.

We should compute the shear stress from the current state.

$$\tau_{\rm slope} = \Delta_{\rho} g D \sin(\alpha) \hat{G},$$

``` {.julia #cat-propagator}
function submarine_transport(input::Input)
  phys_size = (x = input.grid_size[1] * input.phys_scale,
               y = input.grid_size[2] * input.phys_scale)
  box = Transport.Box(input.grid_size, phys_size, input.phys_scale, length(input.facies))

  # This function does not modify the state, rather it transports the
  # sediments in a given product frame and gives a new product frame
  # as output.
  function (s::State, Δ::ProductFrame)  # -> ProductFrame
    output = Array{Float64}(undef, size(Δ.production))
    transport = Transport.transport(input.boundary, box, function (p::Particle)
      z, ∇ = Transport.interpolate(input.boundary, box, s.height, p.position)
      α = atan(abs(∇))
      Ĝ = ∇ / abs(∇)

      τ_wave = input.wave_shear_stress(z)

      Δρ = input.facies[p.facies].excess_density
      D = input.facies[p.facies].grain_size
      g = input.g
      τ_grav = (Δρ * D * g * sin(α)) * Ĝ

      return τ_grav + τ_wave
    end)
    deposit = Transport.deposit(input.boundary, box, output)

    # iterate all particles, transport them and plot
    for (i, (idx, mass)) in enumerate(pairs(Δ.production))
      subgrid_spacing = 1.0 / input.transport_subsample
      subgrid_axis = 0.0:subgrid_spacing:1.0 - subgrid_spacing
      subgrid = product(subgrid_axis, subgrid_axis)
      for (j, dx) in enumerate(subgrid)
        facies_type = idx[1]
        p = (x=(idx[2] + dx[1]) * input.phys_scale, y=(idx[3] + dx[2]) * input.phys_scale)
        θ = input.facies[facies_type].critical_stress
        particle = Particle(p, mass * subgrid_spacing^2, θ, facies_type, nothing) |> transport |> deposit
      end
    end
    return ProductFrame(output)
  end
end

function shear_stress(input::Input)
  # To compute critical shear stress, take weighted average of critical angles
  θ_f = [sin(f.critical_angle) for f in input.facies]
  gradient_stencil = stencil(Float64, Vec2, input.boundary, (3, 3), function (w)
    kernel = [-1 0 1; -2 0 2; -1 0 1] ./ (8 * phys_scale)
    ( x = sum(kernel .* w)
    , y = sum(kernel' .* w) )
  end)

  ∇ = Matrix{Vec2}(undef, input.grid_size...)
  ρD_f = [f.excess_density * f.grain_size for f in input.facies]

  function stress(s::State)
    gradient_stencil(s.height, ∇)
    α = atan.(abs.(∇))

    top_sediment = peek_sediment(s.sediment, input.top_layer_thickness)
    A = sum(top_sediment .* ρD_f .* input.g; dims=3) ./ sum(top_sediment)
    θ = sum(top_sediment .* θ_f; dims=3) ./ sum(top_sediment) 

    gravitational_stress = A .* sin.(α)
    wave_stress = input.wave_shear_stress.(s.height)
    critical_stress = A .* θ
  end
end
```


``` {.julia file=src/CATP.jl}
module CATP

using ..Vectors
using ..Stencil: Periodic
using ..Utility
using ..Stencil
using ..Burgess2013.CA
using ..SedimentStack
using ..Transport

using HDF5

using .Iterators: drop, peel, partition, map, take
<<cat-facies>>
<<cat-input>>
<<cat-state>>
<<cat-frame>>
<<cat-update>>
<<cat-propagator>>
<<ca-prod-propagator>>

end  # CATP
```

## Summary of Carbonate 3D (Warrlich Thesis)
The Carbonate3D model is a (closed-source) competitor to CarboCAT. We may learn from there how to better implement sediment transport.
Three main drivers in this model are the same as in CarboKitten:

- Sea level
- Subsidence
- Sediment Production

### Production
Concerning production, Carbonate3D makes some subtly different choices than were made in CarboCAT. A maximum production rate is found between two depths $z_t$ and $z_b$. From sealevel to $z_t$ there is a linear increase from no production to maximum. This choice could help with numeric stability later on. Below $z_b$ there is an exponential decay. At least in this respect, CarboCAT solves the exponential decay more naturally, using the Bosscher-Schlager model.

Additionally, Warrlich introduces some **stress** factors (he uses the word a bit too much for different concepts), all inducing exponential decay on production rates. For open-water species, these stress factors are summarized into carbonate producing species *liking* to be closer to open water. In effect, using a Gaussian filter to measure the distance to open sea, he stimulates growth at the edges of the platform. No doubt to produce ring-like atol formations later on in the thesis, though it can be argued that this morphology is put in by hand like this.

A second stress factor for open-water species is light intensity reduction by the actual sediment that is being produced, making the water more cloudy. This means that (at least for facies with smaller grain size) sediment production is self-limiting. Also, because sediment is more easily transported away from the platform edge, this makes production more efficient at the edge of the plateau.

A second group of factories is introduced that favour more the shallow conditions on the inside of an atol, although production rates are usually much lower there. The same trick with the Gaussian filter is repeated, this time stimulating growth on the inside of the plateau.

A third factory type: **pelagic production** has much lower production rates but also produces in the open ocean. (pelagic means pertaining to ocean environment)

These three factory types kind of match the three types also found in CarboCAT, even though the model implementation is completely different. What I like about CarboCAT is that it doesn't make the same assumptions, whereas in Carbonate3D a lot of the expected behaviour is more or less programmed in by hand. The attempt of Carbonate3D to include more stress factors (and thus a higher degree of realism) actually decreases the explanatory power of the model.

Carbonate3D adds two more facies for coarse- and fine-grained siliciclastics (floating particles). 

### Dispersal
Sediment dispersal is broken down into:

- Disintegration
- Entrainment (capture of disintegrated sediment and clastics into the water flow)
- Transport
- Redeposition

For each sediment type dispersal is modelled by two parameters: carbonate/siliciclastic ratio and grain/matrix ratio. Here *grain* means coarse particles and *matrix* means fine dust/mud.

Additional sediments are grouped into two classes: clastics from external input (not currently of interest for CarboKitten) and disintegrated older sediments, both parametrized with the same parameters we use for freshly produced facies.

Shallow water platform interior species (algae) are associated with carbonate matrix sediment, same as pelagic production. Products of shallow open water (corals) are mixed size carbonate. Disintegration of older sediments keep these same properties. 

#### Disintegration
Disintegration rate can be a function of depth only, seeing that wave energy is an important factor. Also, it is argued that the shallow areas of the reef interior produce more fine-grained sediment. The second process is driven by gravity due to steep slopes.

#### Entrainment
In this model, sediment is entrained (lifted from the surface and transported) if the shear stress $\tau$ exceeds the critical shear stress $\tau_{e}$.
The shear stress acting on the sediment $\tau_{total}$ can be simplified to a combination of shear stress from currents, waves and slopes:

$$\tau_{\rm total} = \tau_{\rm current} + \tau_{\rm wave} + \tau_{\rm slope}$$

#### Shear stress due to slope

$$\tau_{slope} = \Delta \rho g D \sin( \alpha ) ( \underline{G} )/G$$

$\Delta \rho$ is excess density of the sediment, $g$ is gravitational acceleration, $D$ is grain size and $\alpha$ the max slope angle parallel to the gradient vector $\underline{G}$.

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

