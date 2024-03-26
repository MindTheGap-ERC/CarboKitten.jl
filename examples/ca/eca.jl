# ~/~ begin <<docs/src/stencils.md#examples/ca/eca.jl>>[init]
#| creates: docs/src/fig/eca.png
#| requires: src/Stencil.jl
#| collect: figures

using CarboKitten.BoundaryTrait
using CarboKitten.Stencil
using GLMakie

module ECA
    using CarboKitten.Stencil
    using CairoMakie

    # ~/~ begin <<docs/src/stencils.md#eca-rule>>[init]
    rule(i::Int) = function (foo::AbstractVector{T}) where T <: Integer
        d = foo[1]*4 + foo[2]*2 + foo[3]
        i & (1 << d) == 0 ? 0 : 1
    end
    # ~/~ end

    function eca(r::Int, n::Int, iter::Int)
        y = Array{Int}(undef, n, iter)
        y[:, 1] = zeros(Int, n)
        y[div(n, 2), 1] = 1
        stencil_op = stencil(Int, Periodic{1}, (3,), rule(r))
        for i in 2:iter
            stencil_op(view(y, :, i-1), view(y, :, i))
        end
        y
    end

    function plot()
        fig = Figure(resolution=(1200,400))
        for (idx, r) in enumerate([18, 30, 110])
            ax = Axis(fig[1,idx]; title="rule $(r)", yreversed=true, limits=((1, 256), (1, 128)))
            heatmap!(ax, eca(r, 256, 128); colormap=:Blues)
        end
        save("docs/src/fig/eca.png", fig)
    end
end

ECA.plot()
# ~/~ end