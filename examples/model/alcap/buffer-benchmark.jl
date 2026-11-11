# ~/~ begin <<docs/src/algorithms/circular_stack.md#examples/model/alcap/buffer-benchmark.jl>>[init]
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
# ~/~ end
