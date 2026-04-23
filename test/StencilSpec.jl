# ~/~ begin <<docs/src/stencils.md#test/StencilSpec.jl>>[init]
@testset "CarboKitten.Stencil" begin

using StaticArrays
using CarboKitten
using CarboKitten.Stencil: stencil!

let a = ones(Float64, 10, 10),
    b = zeros(Float64, 10, 10)

    stencil!(Periodic{2}, Size(3, 3), b, a) do a
        sum(a)
    end

    @test all(b .== 9)
end

end
# ~/~ end