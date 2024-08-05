module Denudation

include("./Denudation/Abstract.jl")
include("./Denudation/DissolutionMod.jl")
include("./Denudation/EmpiricalDenudationMod.jl")
include("./Denudation/PhysicalErosionMod.jl")
include("./Denudation/NoDenudationMod.jl")

using .Abstract: DenudationType, denudation, redistribution
using .DissolutionMod: Dissolution
using .EmpiricalDenudationMod: EmpiricalDenudation
using .PhysicalErosionMod: PhysicalErosion
using .NoDenudationMod: NoDenudation

export Dissolution, EmpiricalDenudation, PhysicalErosion, NoDenudation, denudation, redistribution

end
