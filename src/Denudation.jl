module Denudation

include("./Denudation/Abstract.jl")
include("./Denudation/DissolutionMod.jl")
include("./Denudation/EmpericalDenudationMod.jl")
include("./Denudation/PhysicalErosionMod.jl")
include("./Denudation/NoDenudationMod.jl")

using .Abstract: DenudationType, denudation, redistribution
using .DissolutionMod: Dissolution
using .EmpericalDenudationMod: EmpericalDenudation
using .PhysicalErosionMod: PhysicalErosion
using .NoDenudationMod: NoDenudation

export Dissolution, EmpericalDenudation, PhysicalErosion, NoDenudation, denudation, redistribution

end
