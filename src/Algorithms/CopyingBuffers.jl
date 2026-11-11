# ~/~ begin <<docs/src/algorithms/copying_stack.md#src/Algorithms/CopyingBuffers.jl>>[init]
module CopyingBuffers

export push_sediment!, pop_sediment!, peek_sediment

# ~/~ begin <<docs/src/algorithms/copying_stack.md#copying-buffer-impl>>[init]
@inline function shift_down!(col::AbstractMatrix, k::Int)
    s = size(col, 1)
    @inbounds for j in axes(col, 2), i in s:-1:k+1
        col[i, j] = col[i - k, j]
    end
end

@inline function shift_up!(col::AbstractMatrix, k::Int)
    s = size(col, 1)
    @inbounds for j in axes(col, 2), i in 1:s-k
        col[i, j] = col[i + k, j]
    end
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/copying_stack.md#copying-buffer-impl>>[1]
function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
    @assert size(col, 2) == length(parcel) "column $(size(col)) doesn't match parcel $(size(parcel))"
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#push-sediment>>[init]
    mass = sum(parcel)
    if mass == 0.0
        return
    end
    
    if mass > size(col)[1]
        @warn "pushing a very large parcel of sediment: $mass times depositional resolution"
        col .= transpose(parcel ./ mass)
        return
    end
    # ~/~ end
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#push-sediment>>[1]
    bucket = sum(@view col[1, :])
    @assert bucket >= 0.0 && bucket <= 1.0
    # ~/~ end
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#push-sediment>>[2]
    if bucket + mass < 1.0
        @views col[1, :] .+= parcel
        return
    end
    # ~/~ end
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#push-sediment>>[3]
    total = mass
    space = 1.0 - bucket
    @views col[1, :] .+= parcel ./ total .* space
    mass -= space
    n = floor(Int64, mass)
    
    shift_down!(col, n + 1)
    # ~/~ end
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#push-sediment>>[4]
    mass -= n
    for j in axes(col, 2)
        frac = parcel[j] / total
        col[2:n+1, j] .= frac
        col[1, j] = frac * mass
    end
    # ~/~ end
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/copying_stack.md#copying-buffer-impl>>[2]
@inline function pop_fraction!(parcel::AbstractVector{F}, col::AbstractMatrix{F}, mass::F) where F <: Real
    bucket = sum(@view col[1, :])
    if mass == 0 || bucket == 0
        return parcel
    end

    @assert mass < bucket "pop_fraction can only pop from the top cell: $(col), $(mass)"
    r = mass / bucket
    for j in axes(col, 2)
        p = r * col[1, j]
        parcel[j] += p
        col[1, j] -= p
    end
    return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, Δ::F, parcel::AbstractVector{F}) where F <: Real
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#pop-sediment>>[init]
    bucket = sum(@view col[1, :])
    @assert bucket >= 0.0
    
    if Δ < bucket
        parcel .= zero(F)
        return pop_fraction!(parcel, col, Δ)
    end
    # ~/~ end
    # ~/~ begin <<docs/src/algorithms/copying_stack.md#pop-sediment>>[1]
    @views parcel .= col[1, :]
    Δ -= bucket
    n = floor(Int64, Δ)
    
    if n > (size(col)[1] - 2)
        @error "too much material popped of the stack: Δ = $Δ"
        for j in axes(col, 2)
            parcel[j] = sum(@view col[:, j])
        end
        col .= 0.0
        return parcel
    end
    
    for j in axes(col, 2)
        parcel[j] += sum(@view col[2:n+1, j])
    end
    shift_up!(col, n + 1)
    col[end-n:end, :] .= 0
    Δ -= n
    
    return pop_fraction!(parcel, col, Δ)
    # ~/~ end
end

pop_sediment!(col::AbstractMatrix{F}, Δ::F) where F <: Real =
    pop_sediment!(col, Δ, zeros(F, size(col, 2)))
# ~/~ end
# ~/~ begin <<docs/src/algorithms/copying_stack.md#copying-buffer-impl>>[3]
function peek_sediment(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
    if Δ == 0
        return zeros(F, size(col, 2))
    end

    bucket = sum(@view col[1, :])
    if Δ < bucket
        return (Δ / bucket) .* col[1, :]
    end

    parcel = col[1, :]
    Δ -= bucket
    n = min(floor(Int64, Δ), size(col, 1) - 1)

    for j in axes(col, 2)
        parcel[j] += sum(@view col[2:n+1, j])
    end
    Δ -= n

    if n + 2 <= size(col, 1)
        last_bucket = sum(@view col[n+2, :])
        if last_bucket > 0
            @views parcel .+= (Δ / last_bucket) .* col[n+2, :]
        end
    end

    return parcel
end
# ~/~ end
# ~/~ begin <<docs/src/algorithms/copying_stack.md#copying-buffer-impl>>[4]
push_sediment!(buffer::AbstractArray{F, 4}, i::CartesianIndex{2}, parcel::AbstractVector{F}) where F <: Real =
    push_sediment!(view(buffer, :, :, i), parcel)

pop_sediment!(buffer::AbstractArray{F, 4}, i::CartesianIndex{2}, Δ::F, parcel::AbstractVector{F}) where F <: Real =
    pop_sediment!(view(buffer, :, :, i), Δ, parcel)

peek_sediment(buffer::AbstractArray{F, 4}, i::CartesianIndex{2}, Δ::F) where F <: Real =
    peek_sediment(view(buffer, :, :, i), Δ)

function push_sediment!(buffer::AbstractArray{F, 4}, p::AbstractArray{F, 3}) where F <: Real
    @assert size(buffer, 2) == size(p, 1) "shapes for buffer $(size(buffer)) doesn't match p $(size(p))"
    @views for i in CartesianIndices(size(p)[2:3])
        push_sediment!(buffer, i, p[:, i])
    end
end

function pop_sediment!(buffer::AbstractArray{F, 4}, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}) where F <: Real
    @views for i in CartesianIndices(amount)
        pop_sediment!(buffer, i, amount[i], out[:, i])
    end
end

function peek_sediment(buffer::AbstractArray{F, 4}, Δ::F) where F <: Real
    _, f, x, y = size(buffer)
    out = Array{F, 3}(undef, f, x, y)
    for i in CartesianIndices((x, y))
        out[:, i] = peek_sediment(buffer, i, Δ)
    end
    return out
end
# ~/~ end

end # module
# ~/~ end
