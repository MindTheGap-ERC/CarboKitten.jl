# ~/~ begin <<docs/src/components/initial_sediment.md#src/Components/InitialSediment.jl>>[init]
@compose module InitialSediment
    @mixin FaciesBase, WaterDepth, SedimentBuffer
    using Unitful: dimension, NoUnits

    # ~/~ begin <<docs/src/components/initial_sediment.md#initial-sediment>>[init]
    @kwdef struct Facies <: AbstractFacies
        initial_sediment = 0.0u"m"
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/initial_sediment.md#initial-sediment>>[1]
    function initial_sediment(box::Box, facies::AbstractFacies)
        s = facies.initial_sediment
        if s isa AbstractMatrix
            @assert size(s) == box.grid_size
            @assert dimension(eltype(s)) == dimension(u"m")
            return s
        end
        if s isa Quantity
            @assert dimension(s) == dimension(u"m")
            return fill(s, box.grid_size...)
        end

        # s should be callable
        x, y = box_axes(box)
        return s.(x, y')
    end

    function push_initial_sediment!(input::AbstractInput, state::AbstractState)
        s = stack(initial_sediment(input.box, f) for f in input.facies; dims=1)
        push_sediment!(state.sediment_buffer, s ./ input.depositional_resolution .|> NoUnits)
        state.sediment_height .+= sum(s; dims=1)[1,:,:]
    end
    # ~/~ end
end
# ~/~ end
