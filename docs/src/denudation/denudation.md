# Denudation

Denudation can be applied according to three approaches: modelling physical erosion, modelling chemical dissolution and estimating total denudation rates based on chlorine (Cl) isotope data.

1. [Physical Erosion](physical_erosion.md)
2. [Chemical Dissolution](chemical.md)
3. [Empirical Denudation](empirical.md)

## How to use?

In CarboKitten, you choose which type of the three approaches you apply. Please see the examples to do this.

In the `examples/denudation_coupled/` folder of the repository, you find `dissolution_run.jl`, `empirical_run.jl`,`physical_run.jl` for chemical dissolution, empirical denudation or physical denudation, respectively. These files use the sine curve as sea level curve input. The following values of the denudation approach argument are possible:

- `NoDenudation` means no erosion, and is used for debugging only.
- `Dissolution` means chemical dissolution. The default input parameters are: Temperature = 273K, precipitation = 1000 mm/yr, atmospheric $CO_2$ partial pressure = $10^{-1.5}*$ ATM, and reaction rate = 0.002 m/yr.
- `PhysicalErosion` means physical erosion and sediments redistribution. The default parameters is erodability = 0.001 m/yr.
- `EmpiricalDenudation` means total denudation calculated based on empirical relationship by Cl isotope observations. The default input parameter is: precipitation = 1000 mm/yr.

You can change your denudation approach in `const DENUDATION`.

After running the examples, it is recommended to use `denudation_plot.jl` to visualize the sediment accumulation curve and the associated stratigraphic columsn to check the validity. You can change two `const` as you wish.

## Tests for three modes of denudation

In this module, 7 tests are implemented.

Test 1:

```julia
@test sum(denudation_mass_LOW_T) < sum(denudation_mass_HIGH_T)
```

This means that higher temperature would dissolve faster than the colder scenario. It tests the Dissolution mode.

Test2:

```julia
@test sum(denudation_mass_LOW_P) < sum(denudation_mass_HIGH_P)
```

This means more humid scenario has higher denudation rates than the arid scenario. This tests emperical denudation mode.

Test 3:

```julia
@test sum(denudation_mass_phys) > sum(denudation_mass_phys_flat)
```

This means more topography has higher denudation rates than the flatter topography. This tests physical erosion.

Test 4:

```julia
@test sum(denudation_mass_phys) ≈ sum(redistribution_mass)
```

This means in physical erosion mode, the total amount of eroded material = the total amount of redistributed material. In this case, boundary condition of 'Periodic has been used.

Test 5 to 7 are regression tests: the outputs from the module is similar to the values calculated by calculator.

## API

The denudation models all follow the same API.

```{.julia file=src/Denudation/Abstract.jl}
module Abstract

using ...BoundaryTrait: Boundary
using ...Boxes: Box

using Unitful

abstract type DenudationType end

"""
    denudation(box, param, state)


Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.

Returns denudation mass in units of meters.
"""
function denudation(input)

    function (state, water_depth, slope)
        if denudation(input.box, input.denudation, water_depth, slope, input.facies,state) !== nothing
        return denudation(input.box, input.denudation, water_depth, slope, input.facies,state) .* input.time.Δt
        else
        return nothing
        end
    end
end

"""
    denudation(box::Box, param::DenudationType, water_depth, slope, facies)

Computes the amount of denudation. This function is called on a pixel by pixel basis, so all arguments can be assumed to be scalar. The `param` argument should be of a subtype of `DenudationType` containing all the input parameters for this specific denudation model.
"""
function denudation(box::Box, param::DenudationType, water_depth, slope, facies, state)
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

Takes `state`, `water_depth` in meters and `denudation_mass` as a 3D array (facies, x and y coordinates) in units of meters.
"""
function redistribution(input)
    function (state, water_depth, denudation_mass)
        return redistribution(input.box, input.denudation, denudation_mass, water_depth)
    end
end

function redistribution(box::Box, param::DenudationType, denudation_mass, water_depth)
    error("Abstract `redistribution` function called.")
end

end  # module
```

### No Denudation

```{.julia file=src/Denudation/NoDenudationMod.jl}
"""
    module NoDenudation

Doesn't do any denudation: used for testing purposes.
"""
module NoDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Boxes: Box
using Unitful

struct NoDenudation <: DenudationType end

function denudation(box::Box, p::NoDenudation, water_depth::Any, slope, facies, state)
    return nothing
end

function redistribution(box::Box, p::NoDenudation, denudation_mass, water_depth)
    return nothing
end

end
```
