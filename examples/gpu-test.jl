module GPUTest

using oneAPI
using KernelAbstractions
using CarboKitten.BoundaryTrait: Boundary, Periodic, Reflected, Coast
using StaticArrays

@inline modflip(a, l) =
    let b = mod1(a, 2l)
        b > l ? 2l - b + 1 : b
    end

@inline get_bounded(::Type{Periodic{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[mod1.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Reflected{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[modflip.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Coast}, a, i) =
    checkbounds(Bool, a, i) ? a[i] : a[modflip(i[1], size(a)[1]), mod1(i[2], size(a)[2])]

@kernel function life(::Type{BT}, a, b) where {BT}
    i = @index(Global, Cartesian)

    n = 0
    for j in CartesianIndices((3, 3))
        n += get_bounded(BT, a, i + j - CartesianIndex(1, 1))
    end
    n -= a[i]

    b[i] = n == 3 || n == 2 && a[i] == 1
end

@kernel function ca2d5sq(::Type{BT}, in, out, precedence) where {BT}
    i = @index(Global, Cartesian)

    for f in precedence
        n = 0
        for j in CartesianIndices((5, 5))
            n += (get_bounded(BT, in, i + j - CartesianIndex(3, 3)) == f)
        end

        if in[i] == f
            n -= 1
            out[i] = 4 <= n && n <= 10 ? f : 0
            break
        else
            if 6 <= n && n <= 10
                out[i] = f
                break
            end
        end
    end
end

function main()
    dev = CPU()
    a = oneArray(rand(0:3, 256, 256))
    b = oneArray(zeros(Int, 256, 256))
    precedence = collect(1:3)

    backend = get_backend(a)
    for _ in 1:1000
        ca2d5sq(backend, 256)(Periodic{2}, a, b, oneArray(precedence), ndrange=size(a))
        KernelAbstractions.synchronize(backend)
        circshift!(precedence, 1)
        ca2d5sq(backend, 256)(Periodic{2}, b, a, oneArray(precedence), ndrange=size(a))
        KernelAbstractions.synchronize(backend)
        circshift!(precedence, 1)
    end
    a
end

end
