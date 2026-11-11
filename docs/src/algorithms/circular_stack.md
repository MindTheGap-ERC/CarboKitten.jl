# Circular Stack

The [copying stack](copying_stack.md) keeps the head of each column at the first row. Every time the head overflows (or is emptied), all rows below are copied one or more rows down (or up). The circular stack avoids this copying. Instead, it remembers for each column which row is the head. The rows below the head are found by walking downward from it, wrapping around at the end. So, for a buffer of depth $D$, the logical row $r$ (counting from the head) is found at physical row $\textrm{mod1}(h - r + 1, D)$, where $h$ is the head. This means we have something of a [circular buffer](https://en.wikipedia.org/wiki/Circular_buffer).

Pushing sediment that overflows the head moves the head forward, overwriting the deepest rows, which would have been dropped by the copying stack anyway. Popping sediment moves the head back, and zeroes the rows that were popped, which the copying stack would have filled with zeros from below.

The algorithms below follow the copying stack step by step. The result is the same, up to round-off errors and the position of the head.

``` {.julia file=test/Algorithms/CircularBuffersSpec.jl}
module CircularBuffersSpec
using Test
using CarboKitten.Algorithms.CircularBuffers: CircularBuffer, push_sediment!, pop_sediment!, peek_sediment, unroll
using CarboKitten.Algorithms: CopyingBuffers
using Random: Xoshiro

@testset "Algorithms/CircularBuffers" begin
  @testset "column" begin
    col = zeros(Float64, 10, 3)
    head = 1
    parcel = zeros(Float64, 3)
    head = pop_sediment!(col, head, 0.0, parcel)
    @test parcel == [0.0, 0.0, 0.0]
    head = push_sediment!(col, head, [5.0, 0, 0])
    head = pop_sediment!(col, head, 1.5, parcel)
    @test parcel == [1.5, 0.0, 0.0]
    head = push_sediment!(col, head, [0.0, 2.0, 0.0])
    head = pop_sediment!(col, head, 2.0, parcel)
    @test parcel == [0.25, 1.75, 0.0]
    head = pop_sediment!(col, head, 1.5, parcel)
    @test parcel == [1.25, 0.25, 0.0]
    head = pop_sediment!(col, head, 0.0, parcel)
    @test parcel == [0.0, 0.0, 0.0]
  end

  @testset "edge cases" begin
    # pushing more than the buffer can hold fills every row with the parcel's fractions
    col = zeros(Float64, 4, 2)
    head = @test_logs (:warn,) push_sediment!(col, 1, [6.0, 2.0])
    @test all(col .== [0.75 0.25])

    # popping more than the buffer holds returns everything
    col = [0.0 1.0; 0.0 1.0; 1.0 0.0; 1.0 0.0]
    parcel = zeros(Float64, 2)
    @test_logs (:error,) pop_sediment!(col, 4, 10.0, parcel)
    @test parcel == [2.0, 2.0]
    @test all(col .== 0.0)

    # the head wraps around
    col = zeros(Float64, 4, 2)
    head = push_sediment!(col, 4, [1.5, 0.0])
    @test head == 1
    @test col[:, 1] == [0.5, 0.0, 0.0, 1.0]
    @test unroll(col, head) == [0.5 0.0; 1.0 0.0; 0.0 0.0; 0.0 0.0]
  end

  @testset "equivalent to copying stack" begin
    rng = Xoshiro(0)
    for _ in 1:20
      copying = zeros(Float64, 8, 3)
      circular = zeros(Float64, 8, 3)
      head = rand(rng, 1:8)
      for _ in 1:50
        if rand(rng, Bool)
          parcel = rand(rng, 3) .* rand(rng, (0.1, 1.0, 2.0))
          CopyingBuffers.push_sediment!(copying, parcel)
          head = push_sediment!(circular, head, parcel)
        else
          Δ = rand(rng) * rand(rng, (0.1, 1.0, 3.0))
          @test peek_sediment(circular, head, Δ) ≈ CopyingBuffers.peek_sediment(copying, Δ)
          expected = CopyingBuffers.pop_sediment!(copying, Δ)
          parcel = zeros(Float64, 3)
          head = pop_sediment!(circular, head, Δ, parcel)
          @test parcel ≈ expected
        end
        @test unroll(circular, head) ≈ copying
      end
    end
  end

  @testset "grid" begin
    buffer = CircularBuffer{Float64}(10, 3, (5, 5))
    for x in 1:10
      production = rand(3, 5, 5)
      push_sediment!(buffer, production)
    end
    a = peek_sediment(buffer, 1.0)
    @test all(sum(a; dims=1) .≈ 1.0)

    out = zeros(Float64, 3, 5, 5)
    pop_sediment!(buffer, fill(1.0, 5, 5), out)
    @test out ≈ a
  end
end
end
```

## Data structure

The `CircularBuffer` contains the data in the same layout as the copying stack, together with a matrix of heads, one for each column.

``` {.julia #circular-buffer}
struct CircularBuffer{F <: Real}
    data::Array{F, 4}
    heads::Matrix{Int}
end

CircularBuffer{F}(depth::Int, n_facies::Int, grid_size::NTuple{2, Int}) where F <: Real =
    CircularBuffer{F}(zeros(F, depth, n_facies, grid_size...), ones(Int, grid_size...))
```

We need a few helper functions to walk the rows of a column. The `unroll` function returns a column in the layout of the copying stack, which is useful for testing.

``` {.julia #circular-buffer}
@inline row_sum(col::AbstractMatrix, r::Int) = sum(@view col[r, :])
@inline next_row(col::AbstractMatrix, r::Int) = mod1(r + 1, size(col, 1))
@inline prev_row(col::AbstractMatrix, r::Int) = mod1(r - 1, size(col, 1))

unroll(col::AbstractMatrix, head::Int) = col[mod1.(head .- (0:size(col, 1)-1), size(col, 1)), :]
```

## Pushing sediment

The column version of `push_sediment!` takes the column `col`, the current `head` and a `parcel`. It returns the new head. The first steps are identical to those of the copying stack.

``` {.julia #circular-buffer}
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
    <<circular-push-rows>>
end
```

Where the copying stack shifts all rows down by $n + 1$, we move the head forward instead. The $n$ rows we pass are filled with the normalized fractions, the new head with whatever is left.

``` {.julia #circular-push-rows}
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
```

## Popping sediment

Again, `pop_sediment!` returns the new head, writing the popped sediment to `parcel`.

``` {.julia #circular-buffer}
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
    <<circular-pop-rows>>
end
```

Where the copying stack shifts the rows up by $n + 1$, we take the head and the $n$ rows below it, zeroing them and moving the head back as we go. The remainder is taken from the new head.

``` {.julia #circular-pop-rows}
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
```

## Peeking

Peeking walks the same rows as popping, without changing the buffer.

``` {.julia #circular-buffer}
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
```

## Grid

Same as for the copying stack, we have methods that work on a single column of the grid, and methods that loop over all columns.

``` {.julia #circular-buffer}
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
```

## Benchmark

To compare the two implementations, we run the [ALCAP model](../models/alcap.md) once with each buffer. The circular variant of ALCAP only differs in its `Input` and `State` records, which include `CircularSedimentBuffer` instead of the default `CopyingSedimentBuffer`. Both runs should end up with the same sediment.

Most time in a model run is spent on other things than the sediment buffer, like transport. That is why we also time pushing and popping on the final state of each run.

``` {.julia file=examples/model/alcap/buffer-benchmark.jl}
module BufferBenchmark

using Unitful
using CarboKitten
using CarboKitten.Components.Common: Amount
using CarboKitten.Components.SedimentBuffer: push_sediment, pop_sediment
using CarboKitten.SedimentBuffers: peek_sediment
using CarboKitten.Models.ALCAP.Example: FACIES

module CircularALCAP
    using CarboKitten: ALCAP
    using CarboKitten.Components.Common: record_type
    using CarboKitten.Components: CircularSedimentBuffer

    const Input = record_type(:Input, ALCAP, CircularSedimentBuffer)
    const State = record_type(:State, ALCAP, CircularSedimentBuffer)

    initial_state(input) = ALCAP.initial_state(input, State)
    initial_frame(input) = ALCAP.initial_frame(input)
    step!(input) = ALCAP.step!(input)
end

const MODELS = [
    "copying" => (ALCAP, ALCAP.Input),
    "circular" => (CircularALCAP, CircularALCAP.Input)]

make_input(Input, steps) = Input(
    tag="buffer-benchmark",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=steps),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    lithification_time=100.0u"yr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

run(M, input) = run_model((_, _) -> nothing, Model{M}, input)

function best_of(f, n)
    best = (time=Inf, bytes=0)
    for _ in 1:n
        GC.gc()
        t = @timed f()
        if t.time < best.time
            best = (time=t.time, bytes=t.bytes)
        end
    end
    return best
end

function buffer_ops(input, state)
    push! = push_sediment(input)
    pop! = pop_sediment(input)
    deposit = fill(0.4u"m", length(input.facies), input.box.grid_size...)
    amount = fill(1.2u"m", input.box.grid_size...)
    out = zeros(Amount, length(input.facies), input.box.grid_size...)
    return () -> for _ in 1:100
        push!(state, deposit)
        pop!(state, amount, out)
    end
end

function main(; steps=1000, repeats=3)
    states = Dict()
    println("ALCAP, 100×50 grid, $(steps) steps, best of $(repeats)")
    for (name, (M, Input)) in MODELS
        run(M, make_input(Input, 10))  # compile
        input = make_input(Input, steps)
        model = best_of(() -> states[name] = run(M, input), repeats)
        ops = best_of(buffer_ops(input, deepcopy(states[name])), repeats)
        println(rpad(name, 10),
            "model run: ", lpad(round(model.time; digits=2), 7), " s  ",
            lpad(Base.format_bytes(model.bytes), 10), "    ",
            "100 × push+pop: ", lpad(round(ops.time * 1000; digits=1), 7), " ms  ",
            lpad(Base.format_bytes(ops.bytes), 10))
    end

    copying, circular = states["copying"], states["circular"]
    same_thickness = copying.sediment_thickness ≈ circular.sediment_thickness
    same_buffer = peek_sediment(copying, 50.0) ≈ peek_sediment(circular, 50.0)
    println("same result: thickness $(same_thickness), buffer contents $(same_buffer)")
end

end

BufferBenchmark.main()
```

```@raw html
<details><summary>CircularBuffers module</summary>
```

``` {.julia file=src/Algorithms/CircularBuffers.jl}
module CircularBuffers

export CircularBuffer, push_sediment!, pop_sediment!, peek_sediment

<<circular-buffer>>

end # module
```

```@raw html
</details>
```
