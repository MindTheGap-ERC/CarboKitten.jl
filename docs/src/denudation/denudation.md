# Denudation

Denudation could be achieved by three ways: modelling physical erosion, modelling chemical dissolution and estimating total denudation rates based on chlorine (Cl) isotope data.

1. [Physical Erosion](physical_erosion.md)
2. [Chemical Dissolution](chemical.md)
3. [Empirical Denudation](empirical.md)

## How to use?
In CarboKitten, you could choose which type of the three you would like to attempt. To do this you could simply change the `erosion_type` in the input.

Example: in examples, you find `caps_miller_diss.jl`, `caps_miller_emp.jl`, `caps_miller_phys.jl`, for chemical dissolution, empirical denudation or physical denudation, respectively. This file uses the [miller_phanerozoic_2005] (@cite) cure as sea level curve input. You could try different erosion types by changing the `erosion_type`:

- `NoDenudation` means no erosion, and is used for debugging only.
- `Dissolution` means chemical dissolution
- `PhysicalErosion` means physical erosion and sediments redistribution
- `EmpericalDenudation` means total denudation calculated based on emperical relationship by Cl isotope observations.

## API

The different denudation models all follow the same API.

``` {.julia file=src/Denudation/Abstract.jl}
module Abstract

using ...BoundaryTrait: Boundary

abstract type DenudationType end

"""
    denudation(box, param, state)

FIXME Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(::Type{BT}, param::DenudationType, water_depth, slope, facies) where {BT<:Boundary}
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

FIXME
"""
function redistribution(::Type{BT}, param::DenudationType, water_depth, slope, facies) where {BT<:Boundary}
    error("Abstract `redistribution` function called.")
end

end  # module
```

### No Denudation

``` {.julia file=src/Denudation/NoDenudationMod.jl}
"""
    module NoDenudation

Doesn't do any denudation: used for testing purposes.
"""
module NoDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Config: Box
using Unitful

struct NoDenudation <: DenudationType end

function denudation(box::Box, p::NoDenudation, water_depth::Any, slope, facies)
    return nothing
end

function redistribution(box::Box, p::NoDenudation, water_depth, slope, facies)
    return nothing
end

end
```
