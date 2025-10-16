# ~/~ begin <<docs/src/architecture.md#src/RunModel.jl>>[init]
module RunModel

import ..CarboKitten: n_steps, run_model, get_logger, Model, AbstractInput, AbstractState
using ProgressLogging
using Logging: with_logger

# ~/~ begin <<docs/src/architecture.md#run-model>>[init]
"""
    run_model(f, ::Type{Model{M}}, input::AbstractInput) where M
    run_model(f, ::Type{Model{M}}, input::AbstractInput, state::AbstractState) where M

Run a model and send `Frame`s to callback `f`. The type parameter `M` should be a model,
being a module with both `initial_state!` and `step!` defined.

The second version with explicit state is used by `H5Writer` so we can perform an additional
action between creating the initial state and starting the model run (saving metadata to the
HDF5 file).
"""
run_model(f, ::Type{Model{M}}, input::AbstractInput) where M =
    run_model(f, Model{M}, input, M.initial_state(input))

function run_model(f, ::Type{Model{M}}, input::AbstractInput, state::AbstractState) where M
    logger = get_logger(input)
    with_logger(logger) do
        step! = M.step!(input)

        @progress for w = 1:n_steps(input)
            f(w, step!(state))
        end
    end
end
# ~/~ end

end
# ~/~ end
