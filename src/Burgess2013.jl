# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013.jl>>[init]
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
# ~/~ end