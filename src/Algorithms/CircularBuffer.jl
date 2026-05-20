# ~/~ begin <<docs/src/algorithms/circular_stack.md#src/Algorithms/CircularBuffer.jl>>[init]
module CircularBuffer

using KernelAbstractions
using StaticArrays
using KernelAbstractions.Extras.LoopInfo: @unroll
using Adapt

# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[init]
struct CircularBuffer{SB,BH}
    data::SB
    heads::BH
end

@adapt_structure CircularBuffer

function initial_state(input)
    n_f = length(input.facies)
    n_b = input.sediment_buffer_size
    n_g = input.box.grid_size
    SedimentBuffer(
        zeros(Float32, n_b, n_f, n_g...),
        zeros(Int32, n_g...))
end
# ~/~ end

@inline function pop_sediment!(
    ::Size{Sz},
    buffer::AbstractArray{F, 1},
    stack::AbstractArray{F,2},
    head::I,
    mass::F) where {Sz, F <: Real, I <: Integer}

    (D, N) = Sz

    exported = MVector{N, F}((0.0 for _ in 1:N)...)

    while mass > 0.0
        bucket::F = zero(F)
        for f in 1:N
            bucket += stack[head, f, i, j]
        end
        if bucket > 0.0
            x = min(bucket, mass)
            for f in 1:N
                xf = x * stack[head, f] / bucket
                stack[head, f] -= xf
                exported[f] += xf
            end
            mass -= x
        else
            mass = 0.0
        end
        head = mod1(head - 1, D)
    end

    for f in 1:N
        buffer[f] = exported[f]
    end

    return mod1(head + 1, D)
end


@kernel function pop_sediment_k(
    ::Size{Sz},
    buffer::AbstractArray{F, 3},
    stack::AbstractArray{F, 4},
    heads::AbstractArray{I, 2},
    @Const mass::AbstractArray{F, 2}) where {Sz, F <: Real, I <: Integer}

    (D, N) = Sz
    i, j = @index(Global, NTuple)

    m::F = mass[i, j]
    # There is a macro `@MVector zeros(N)` but it has a scoping issue
    # not recognizing symbol N
    exported::MVector{N, F} = MVector{N, F}((0.0 for _ in 1:N)...)
    head::I = heads[i, j]

    while m > 0.0
        bucket::F = zero(F)
        @unroll for f in 1:N
            bucket += stack[head, f, i, j]
        end
        if bucket > 0.0
            x = min(bucket, m)
            @unroll for f in 1:N
                xf = x * stack[head, f, i, j] / bucket
                stack[head, f, i, j] -= xf
                exported[f] += xf
            end
            m = m - x
        else
            m = 0.0
        end
        head = mod1(head - 1, D)
    end

    @unroll for f in 1:N
        buffer[f, i, j] += exported[f]
    end

    heads[i, j] = mod1(head + 1, D)
end


@inline function push_sediment!(
    ::Size{Sz},
    stack::AbstractArray{F, 2},
    head::I,
    sediment::AbstractArray{F, 1}) where {Sz, F <: Real, I <: Integer}

    (D, N) = Sz
    mass::F = +((sediment[f, i, j] for f in 1:N)...)
    fractions = SVector{3,F}((sediment[f] / mass for f in 1:N)...)
    while mass > 0.0
        bucket::F = zero(F)
        for f in 1:N
            bucket += stack[head, f]
        end
        x::F = min(1.0 - bucket, mass)
        for f in 1:N
            stack[head, f] += x * fractions[f]
        end
        mass = mass - x
        head = mod1(head + 1, D)
        for f in 1:N
            stack[head, f] = 0.0
        end
    end
    return mod1(head - 1, D)
end


@kernel function push_sediment_k(
    sz::Size{Sz},
    stack::AbstractArray{F, 4},
    heads::AbstractArray{I,2},
    @Const sediment::AbstractArray{F, 3}) where {Sz, F, I}

    (D, N) = Sz
    i, j = @index(Global, NTuple)
    # heads[i, j] = push_sediment!(sz, @view(stack[:, :, i, j]), heads[i, j], @view(sediment[:, i, j]))
    head::Int32 = heads[i, j]

    # mass::F = zero(F)
    # @unroll for f in 1:N
    #     mass += sediment[f, i, j]
    # end
    # # CUDA: Somehow, the above code triggers an illegal memory access error.
    # # The original had `mass::F = sum(@view sediment[:, i, j])` but this
    # # uses dynamic code to determine the length of the loop, while we have
    # # that information statically available.

    mass::F = +((sediment[f, i, j] for f in 1:N)...)
    fractions = SVector{N, F}((sediment[f, i, j] / mass for f in 1:N)...)
    while mass > 0.0
        # bucket::F = +((stack[head, f, i, j] for f in 1:N)...)
        # # Somehow, the above code triggers an LLVM error
        bucket::F = zero(F)
        @unroll for f in 1:N
            bucket += stack[head, f, i, j]
        end

        x = min(1.0 - bucket, mass)
        @unroll for f in 1:N
            stack[head, f, i, j] += x * fractions[f]
        end
        mass = mass - x
        head = mod1(head + 1, D)
        @unroll for f in 1:N
            stack[head, f, i, j] = 0.0
        end
    end
    heads[i, j] = mod1(head - 1, D)
end

end  # module
# ~/~ end
