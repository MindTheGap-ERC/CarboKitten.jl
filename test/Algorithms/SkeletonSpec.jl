# ~/~ begin <<docs/src/algorithms/skeleton.md#test/Algorithms/SkeletonSpec.jl>>[init]
module SkeletonSpec
using Test
using CarboKitten.Algorithms: skeleton

centered_grid(n) = ((i - (n + 1) / 2, j - (n + 1) / 2) for i in 1:n, j in 1:n)
disc(n, r) = [hypot(x, y) <= r for (x, y) in centered_grid(n)]
annulus(n, r1, r2) = [r1 <= hypot(x, y) <= r2 for (x, y) in centered_grid(n)]

function n_components(n_vertices, edges)
    parent = collect(1:n_vertices)
    root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
    for (a, b) in eachrow(edges)
        parent[root(a)] = root(b)
    end
    return count(i -> root(i) == i, 1:n_vertices)
end

"""
Returns the number of connected components and the number of independent
loops in the skeleton of `bitmap`.
"""
function topology(bitmap; minwidth=1)
    vertices, edges = skeleton(bitmap; minwidth=minwidth)
    v = length(vertices)
    e = size(edges, 1)
    c = n_components(v, edges)
    return (components=c, loops=e - v + c)
end

@testset "Algorithms/Skeleton" begin
    @testset "graph structure" begin
        vertices, edges = skeleton(annulus(41, 8, 18); minwidth=1)
        @test all(1 .<= edges .<= length(vertices))
        # edges only connect vertices in adjacent rows
        @test all(vertices[b][1] - vertices[a][1] == 1 for (a, b) in eachrow(edges))
    end

    @testset "topology" begin
        @test topology(disc(41, 15)) == (components=1, loops=0)
        @test topology(annulus(41, 8, 18)) == (components=1, loops=1)
        @test topology(annulus(41, 15, 17)) == (components=1, loops=1)
        @test topology(hcat(disc(21, 8), disc(21, 8))) == (components=2, loops=0)
    end

    @testset "minwidth" begin
        # The walls of the annulus are about 11 cells wide. A larger `minwidth`
        # removes the sides, leaving the top and bottom caps as two trees.
        @test topology(annulus(41, 8, 18); minwidth=12) == (components=2, loops=0)
        @test isempty(skeleton(annulus(41, 8, 18); minwidth=40)[1])
    end

    @testset "empty input" begin
        vertices, edges = skeleton(falses(5, 5); minwidth=1)
        @test isempty(vertices)
        @test size(edges) == (0, 2)
    end
end
end
# ~/~ end
