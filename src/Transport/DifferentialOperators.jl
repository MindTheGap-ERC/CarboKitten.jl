# ~/~ begin <<docs/src/finite-difference-transport.md#src/Transport/DifferentialOperators.jl>>[init]
module DifferentialOperators

# ~/~ begin <<docs/src/finite-difference-transport.md#differential-operators>>[init]
central_difference(::Type{Val{1}}, a::AbstractMatrix, dx) =
    (a[3, 2] - a[1, 2]) / (2dx)

central_difference(::Type{Val{2}}, a::AbstractMatrix, dx) =
    (a[2, 3] - a[2, 1]) / (2dx)
# ~/~ end
# ~/~ begin <<docs/src/finite-difference-transport.md#differential-operators>>[1]
upwind(v::T, a1, a2, a3, dx) where {T} =
    if v < zero(T)
        v * (a3 - a2) / dx
    else
        v * (a2 - a1) / dx
    end

upwind(::Type{Val{1}}, v, a::AbstractMatrix, dx) =
    upwind(v, a[:, 2]..., dx)

upwind(::Type{Val{2}}, v, a::AbstractMatrix, dx) =
    upwind(v, a[2, :]..., dx)
# ~/~ end
# ~/~ begin <<docs/src/finite-difference-transport.md#differential-operators>>[2]
laplacian(a::AbstractMatrix, dx) =
    (a[1, 2] + a[2, 1] + a[3, 2] + a[2, 3] - 4 * a[2, 2]) / dx^2
# ~/~ end

end
# ~/~ end