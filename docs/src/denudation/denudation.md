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
