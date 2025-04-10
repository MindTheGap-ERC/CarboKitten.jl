# Running a Model
Running a model, we store the changes in each time step. The `Frame` type stores the production, disintegration and deposited material after transport.

``` {.julia file=src/Components/Models.jl}
@compose module Models
    using ..Common
    using ..FaciesBase: n_facies
    using ..TimeIntegration: n_writes
    import ...CarboKitten: run_model

    using Unitful
    using ProgressLogging

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

    Run a model and send `Frame`s to callback `f`.
    """
    run_model(f, ::Type{Model{M}}, input::AbstractInput) where M =
        run_model(f, Model{M}, input, M.initial_state(input))

    function run_model(f, ::Type{Model{M}}, input::AbstractInput, state::AbstractState) where M
        step! = M.step!(input)

        @progress for w = 1:n_writes(input)
            df = zeros(Frame, input)
            for n = 1:input.time.write_interval
                increment!(df, step!(state))
            end
            f(w, df)
        end
    end

end
```

