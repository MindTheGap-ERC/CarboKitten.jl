# ~/~ begin <<docs/src/components/run_model.md#src/Components/Models.jl>>[init]
@compose module Models
    @mixin FaciesBase, TimeIntegration
    using ..Common
    using ..FaciesBase: n_facies
    using ..TimeIntegration: n_writes

    import ...CarboKitten: run_model

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

end
# ~/~ end
