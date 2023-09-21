# About
CarboCAT is primarily based on a very simple cellular automaton (CA). We may explore this CA as a first step in implementing the model in Julia.

# Overview

The CarboCAT model [B13: @Burgess2013] consists of several components, many of which are optional or contain optional levels of complexity.

1. Species **habitation**: an algorithm is in place to evolve the locality of a number of factory species.
2. Sediment **production**: each species will produce sediment according to some model.
3. Transport: sediment may be **transported** from a production site to elsewhere due to gravity, waves or other types of mixing.
4. Erosion: sediment may **erode** depending on local circumstances or sediment type.
5. Compactification: different types of sediment may respond to **compression** forces differently.

These processes describe the intrinsic properties of the model. Any parameters that change these processes will be referred to as **model parameters**. Next to that, there are some extrinsic parameters that change the specific output of a model: the initial depth of the sea bed (also known as *bathymetry*) and variation in sea level (including subsidence). These we call **input parameters**.

## Carbonate production

By itself, a sediment production model is enough to model a cross-section of a carbonate platform [BS92: @Bosscher1992]. As a first step, we have [reproduced some results of BS92](bosscher-1992.md). Using a reasonably simple approximation of a growth rate as

[$$\partial_t h = - g_m \tanh \left[\frac{I_0}{I_k} \exp(-k * (h - s(t))\right],$$]{#eq:growth-rate-eqn}

where $h$ is the depth of the sea floor, $g_m$ is the maximum growth rate, $I_0$ the surface light intensity, $I_k$ the saturating light intensity, $k$ the extinction coefficient, and $s$ the (extrinsic) sea-level. In one example given by BS92, we arrived at the following profile.

!include docs/fig/bs92-fig8.html

## Species habitation

These species can be anything, just remember that they are the progenitor of some (limestone) facies type. In the original 2013 model, this stage is implemented by a celullar automaton (or CA). The CA has the nice property of giving pseudo-random output with at least some degree of coherence. There is no physical basis to the CA model, but neither is there very much data to test a physical model against.

We have [implemented the CA used in Burgess 2013](carbocat-ca.md). Using three species with identical 4-6-10-10 rules (survival between 4 to 10 neighbours, birth between 6-10 live neighbours).

![](fig/burgess2013-fig3.svg)

An interesting question is under what rules is this CA stable (i.e. keeps evolving)?

## Combination

The minimal Carbocat model would consist of only species habitation and production.



:::details
### Some submodules

``` {.julia file=src/Burgess2013.jl}
module Burgess2013

include("Burgess2013/Config.jl")
include("Burgess2013/CA.jl")
include("Burgess2013/Production.jl")
include("Burgess2013/Transport.jl")

using .CA
using .Config
using .Production

export production_rate, run_ca, Facies

end
```

``` {.julia #ck-types}
export Product

struct Product
    species::Int
    amount::Float64
end

Base.zero(::Type{Product}) = Product(0, 0.0)
```

``` {.julia file=src/Burgess2013/Config.jl}
module Config

export Facies, MODEL1

struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

MODEL1 = [
    Facies((4, 10), (6, 10), 500, 0.8, 300),
    Facies((4, 10), (6, 10), 400, 0.1, 300),
    Facies((4, 10), (6, 10), 100, 0.005, 300)
]

end
```
:::


``` {.julia file=src/Burgess2013/Production.jl}
module Production

export production_rate

using ..Config: Facies

<<carbonate-production>>

function production_rate(insolation::Float64, facies::Facies, water_depth::Float64)
    gₘ = facies.maximum_growth_rate
    I₀ = insolation
    Iₖ = facies.saturation_intensity
    w = water_depth
    k = facies.extinction_coefficient
    return w > 0.0 ? gₘ * tanh(I₀/Iₖ * exp(-w * k)) : 0.0
end

end
```

### Crowding
In crowded areas carbonate production rates are reduced. For cells where

$$n_{min} \le n \le n_{opt}$$

and

$$n_{opt} \le n \le n_{max}$$

($n_{min}$ and $n_{max}$ for living cells are 4 and 10)

we have a linear increase and linear decrease of production rate (i.e. a triangle function).

## Subsidence
Subsidence refers to gradual lowering or lifting of the underlying floor bed. This could be either sea level rise due to climate change, but also tectonic activity. Sea level also changes according to a given recipe with say three sinusoidals (e.g. Milankovich cycles). When a cell gets "subaerial exposure", i.e. the water level drops below the cell elevation (stupid jargon), accumulation stops and the cell becomes dormant. On reflooding, the cell resumes activity. From the text it is not entirely clear, but I think deactivated cells don't take part in the CA, so they count as dead neighbours.

- [reference on accommodation](http://strata.uga.edu/sequence/accommodation.html), also links to a model called [SedFlux](https://github.com/mcflugen/sedflux).

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.


## Steps

1. Update waterdepth, given subsidence
2. Update sea level elevation, given eustatics
3. Run CA
4. Compute thickness of carbonate production
5. Compute sediment transport

We need to keep a state with the following components:

- height map
- species
- global time, implying:
  - sea level and subsidence

Every cycle we may export a layer of sediment to archive.

# References
