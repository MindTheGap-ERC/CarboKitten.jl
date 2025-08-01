# ~/~ begin <<docs/src/components/run_model.md#src/Components/Models.jl>>[init]
@compose module Models
    @mixin FaciesBase, TimeIntegration
    using ..Common
    using ..FaciesBase: n_facies
    using ..TimeIntegration: n_writes

    import ...CarboKitten: run_model
    using ...OutputData

    using Unitful
    using ProgressLogging

    Base.zeros(::Type{Frame}, input::AbstractInput) = Frame(
        disintegration=zeros(Sediment, n_facies(input), input.box.grid_size...),
        production=zeros(Sediment, n_facies(input), input.box.grid_size...),
        deposition=zeros(Sediment, n_facies(input), input.box.grid_size...))

    function increment!(a::Frame, b::Frame)
        if !isnothing(b.disintegration)
            a.disintegration .+= b.disintegration
        end
        if !isnothing(b.production)
            a.production .+= b.production
        end
        if !isnothing(b.deposition)
            a.deposition .+= b.deposition
        end
    end

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
        step! = M.step!(input)

        @progress for w = 1:n_writes(input)
            df = zeros(Frame, input)
            increment!(df, step!(state))
            f(w, df)
        end
    end

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
end
# ~/~ end
