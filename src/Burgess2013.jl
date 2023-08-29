# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013.jl>>[init]
module Burgess2013

module Types
    # ~/~ begin <<docs/src/carbocat-transport.md#ck-types>>[init]
    struct Deposit{N}
        amount::NTuple{N, Float64}
    end

    Base.zero(::Type{Deposit{N}}) where {N} =
        Deposit{N}(ntuple(_ -> zero(Float64), N))
    # ~/~ end
    # ~/~ begin <<docs/src/carbocat.md#ck-types>>[0]
    export Product

    struct Product
        species::Int
        amount::Float64
    end

    Base.zero(::Type{Product}) = Product(0, 0.0)
    # ~/~ end
end

include("Burgess2013/Config.jl")
include("Burgess2013/CA.jl")
include("Burgess2013/Production.jl")
include("Burgess2013/Transport.jl")

# ~/~ begin <<docs/src/carbocat-cpt.md#carbocat-composite>>[init]
function run()
    # Run the CA for 10 generations as a warm-up
    species = Iterators.drop(CA.run(Reflected{2}, rand(0:3, 50, 50), 3), 10)
    # Initial depth runs from 0 to 150
    height = repeat(collect(0:49) .* 3.0, 1, 50)
    sealevel(t) = 0.0
    Δt = 1000.0

    for (time_index, gen) in enumerate(species)
        t = time_index * Δt
        production = Δt .* Production.production_rate(2000.0, Config.model1, height .- sealevel(t))

    end
end
# ~/~ end

end
# ~/~ end