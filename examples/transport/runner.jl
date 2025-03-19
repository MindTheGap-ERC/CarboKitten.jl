# ~/~ begin <<docs/src/finite-difference-transport.md#examples/transport/runner.jl>>[init]
module Runner
using CarboKitten
using ProgressLogging

n_steps(input) = input.time.steps

function run_model(f, ::Type{Model{M}}, input) where {M}
    state = M.initial_state(input)
    f(0, state)

    @progress for i = 1:n_steps(input)
        M.step!(input, state)
        f(i, state)
    end

    return state
end

do_nothing(_i, _s) = ()

run_model(::Type{Model{M}}, input) where {M} = run_model(do_nothing, Model{M}, input)
end
# ~/~ end