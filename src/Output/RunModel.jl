# ~/~ begin <<docs/src/output/abstract.md#src/Output/RunModel.jl>>[init]
module RunModel

import ...CarboKitten: run_model, Model
using ..Abstract 

# ~/~ begin <<docs/src/output/abstract.md#run-model-output>>[init]
"""
    run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where M

Run a model and save the output to `output`.
"""
function run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where M
    state = M.initial_state(input)
    write_state = state_writer(input, output)
    write_frame = frame_writer(input, output)

    # create a group for every output item
    for (k, v) in input.output
        add_data_set(output, k, v)
    end
    write_state(1, state)

    run_model(Model{M}, input, state) do w, df
        # write_frame chooses to advance in a dataset
        # or just to increment on the current frame
        write_frame(w, df)
        # write_state only writes one in every write_interval
        # and does no accumulation
        write_state(w+1, state)
    end

    return output
end
# ~/~ end

end
# ~/~ end