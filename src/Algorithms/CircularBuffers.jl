# ~/~ begin <<docs/src/algorithms/circular_stack.md#src/Algorithms/CircularBuffers.jl>>[init]
module CircularBuffers

export CircularBuffer, push_sediment!, pop_sediment!, peek_sediment

# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[init]
struct CircularBuffer{F <: Real}
    data::Array{F, 4}
    heads::Matrix{Int}
end

CircularBuffer{F}(depth::Int, n_facies::Int, grid_size::NTuple{2, Int}) where F <: Real =
    CircularBuffer{F}(zeros(F, depth, n_facies, grid_size...), ones(Int, grid_size...))
# ~/~ end
# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[1]
@inline row_sum(col::AbstractMatrix, r::Int) = sum(@view col[r, :])
@inline next_row(col::AbstractMatrix, r::Int) = mod1(r + 1, size(col, 1))
@inline prev_row(col::AbstractMatrix, r::Int) = mod1(r - 1, size(col, 1))

unroll(col::AbstractMatrix, head::Int) = col[mod1.(head .- (0:size(col, 1)-1), size(col, 1)), :]
# ~/~ end
# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[2]
function push_sediment!(col::AbstractMatrix{F}, head::Int, parcel::AbstractVector{F}) where F <: Real
    @assert size(col, 2) == length(parcel) "column $(size(col)) doesn't match parcel $(size(parcel))"
    mass = sum(parcel)
    if mass == 0.0
        return head
    end

    if mass > size(col, 1)
        @warn "pushing a very large parcel of sediment: $mass times depositional resolution"
        col .= transpose(parcel ./ mass)
        return head
    end

    bucket = row_sum(col, head)
    @assert bucket >= 0.0 && bucket <= 1.0

    if bucket + mass < 1.0
        @views col[head, :] .+= parcel
        return head
    end

    total = mass
    space = 1.0 - bucket
    @views col[head, :] .+= parcel ./ total .* space
    mass -= space
    n = floor(Int64, mass)
    # ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-push-rows>>[init]
    for _ in 1:n
        head = next_row(col, head)
        for j in axes(col, 2)
            col[head, j] = parcel[j] / total
        end
    end
    
    head = next_row(col, head)
    mass -= n
    for j in axes(col, 2)
        col[head, j] = parcel[j] / total * mass
    end
    return head
    # ~/~ end
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[3]
@inline function pop_fraction!(parcel::AbstractVector{F}, col::AbstractMatrix{F}, head::Int, mass::F) where F <: Real
    bucket = row_sum(col, head)
    if mass == 0 || bucket == 0
        return parcel
    end

    @assert mass < bucket "pop_fraction can only pop from the top cell: $(col), $(mass)"
    r = mass / bucket
    for j in axes(col, 2)
        p = r * col[head, j]
        parcel[j] += p
        col[head, j] -= p
    end
    return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, head::Int, Δ::F, parcel::AbstractVector{F}) where F <: Real
    bucket = row_sum(col, head)
    @assert bucket >= 0.0

    parcel .= zero(F)
    if Δ < bucket
        pop_fraction!(parcel, col, head, Δ)
        return head
    end

    Δ -= bucket
    n = floor(Int64, Δ)

    if n > size(col, 1) - 2
        @error "too much material popped of the stack: Δ = $Δ"
        for j in axes(col, 2)
            parcel[j] = sum(@view col[:, j])
        end
        col .= 0.0
        return head
    end
    # ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-pop-rows>>[init]
    for _ in 0:n
        for j in axes(col, 2)
            parcel[j] += col[head, j]
            col[head, j] = zero(F)
        end
        head = prev_row(col, head)
    end
    Δ -= n
    
    pop_fraction!(parcel, col, head, Δ)
    return head
    # ~/~ end
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[4]
function peek_sediment(col::AbstractMatrix{F}, head::Int, Δ::F) where F <: Real  # -> Vector{F}
    parcel = zeros(F, size(col, 2))
    if Δ == 0
        return parcel
    end

    bucket = row_sum(col, head)
    if Δ < bucket
        @views parcel .= (Δ / bucket) .* col[head, :]
        return parcel
    end

    Δ -= bucket
    n = min(floor(Int64, Δ), size(col, 1) - 1)
    for _ in 0:n
        @views parcel .+= col[head, :]
        head = prev_row(col, head)
    end
    Δ -= n

    if n + 2 <= size(col, 1)
        last_bucket = row_sum(col, head)
        if last_bucket > 0
            @views parcel .+= (Δ / last_bucket) .* col[head, :]
        end
    end

    return parcel
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/circular_stack.md#circular-buffer>>[5]
function push_sediment!(buffer::CircularBuffer{F}, i::CartesianIndex{2}, parcel::AbstractVector{F}) where F <: Real
    buffer.heads[i] = push_sediment!(view(buffer.data, :, :, i), buffer.heads[i], parcel)
end

function pop_sediment!(buffer::CircularBuffer{F}, i::CartesianIndex{2}, Δ::F, parcel::AbstractVector{F}) where F <: Real
    buffer.heads[i] = pop_sediment!(view(buffer.data, :, :, i), buffer.heads[i], Δ, parcel)
    return parcel
end

peek_sediment(buffer::CircularBuffer{F}, i::CartesianIndex{2}, Δ::F) where F <: Real =
    peek_sediment(view(buffer.data, :, :, i), buffer.heads[i], Δ)

function push_sediment!(buffer::CircularBuffer{F}, p::AbstractArray{F, 3}) where F <: Real
    @assert size(buffer.data, 2) == size(p, 1) "shapes for buffer $(size(buffer.data)) doesn't match p $(size(p))"
    @views for i in CartesianIndices(buffer.heads)
        push_sediment!(buffer, i, p[:, i])
    end
end

function pop_sediment!(buffer::CircularBuffer{F}, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}) where F <: Real
    @views for i in CartesianIndices(amount)
        pop_sediment!(buffer, i, amount[i], out[:, i])
    end
end

function peek_sediment(buffer::CircularBuffer{F}, Δ::F) where F <: Real
    _, f, x, y = size(buffer.data)
    out = Array{F, 3}(undef, f, x, y)
    for i in CartesianIndices((x, y))
        out[:, i] = peek_sediment(buffer, i, Δ)
    end
    return out
end
# ~/~ end

end # module
# ~/~ end
