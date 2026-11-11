# ~/~ begin <<docs/src/components/initial_sediment.md#src/Components/InitialSediment.jl>>[init]
@component module InitialSediment
    @mixin FaciesBase, WaterDepth, SedimentBuffer
    using ..Common
    import ..Common: write_header
    using ..SedimentBuffer: push_sediment
    using Unitful: dimension, NoUnits

    # InitialSediment has no input parameters of its own, but needs a part to
    # dispatch `write_header` on.
    struct Input <: InputPart end

    # ~/~ begin <<docs/src/components/initial_sediment.md#initial-sediment>>[init]
    @kwdef struct Facies <: FaciesPart
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
        push_sediment(input)(state, s)
    end
    
    Common.initialize!(::Input, input::AbstractInput, state::AbstractState) =
        push_initial_sediment!(input, state)
    # ~/~ end

    function write_header(::Input, input::AbstractInput, output::AbstractOutput)
        for (i,f) in enumerate(input.facies)
            set_attribute(output, "facies$(i)/initial_sediment", initial_sediment(input.box, f) |> in_units_of(u"m"))
        end
    end
end
# ~/~ end
