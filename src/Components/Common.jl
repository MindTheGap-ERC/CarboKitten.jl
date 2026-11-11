# ~/~ begin <<docs/src/components/components.md#src/Components/Common.jl>>[init]
module Common
export @u_str, Quantity, Amount, Time, Location, Rate, Intensity, Height, Sediment
export AbstractFacies, AbstractInput, AbstractState, AbstractOutput, set_attribute
export InputPart, StatePart, FaciesPart, write_header
export OutputSpec
export Box, box_axes, Boundary, Coast, Periodic, Reflected, TimeProperties
export in_units_of
export Model
export record_type
export Size
export n_steps
export AbstractProduction, NoProduction, BenthicProduction, PelagicProduction, production_profile, InterpolatedProduction

using ModuleMixins.Procedural: record_type, parts, map_components
using Unitful
using StaticArrays

using ...CarboKitten: Model, AbstractInput, AbstractState, AbstractFacies, AbstractOutput, OutputSpec, set_attribute, n_steps
using ...CarboKitten: InputPart, StatePart, FaciesPart
using ...BoundaryTrait
using ...Config: TimeProperties
using ...Boxes: Box, box_axes
using ...Utility: in_units_of
using ...Production

const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Height = typeof(1.0u"m")
const Location = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")
const Intensity = typeof(1.0u"W/m^2")
const Sediment = typeof(1.0u"m")

"""
    write_header(part, input, output)

Write the attributes of a single component to the output header. Components
extend this method for their own `Input` part.
"""
write_header(part, input::AbstractInput, output::AbstractOutput) = nothing

"""
    write_header(input, output)

Write the header for all components of a model.
"""
write_header(input::AbstractInput, output::AbstractOutput) =
    foreach(p -> write_header(p, input, output), parts(input))

"""
    initial_state(part, input)

The initial `State` part of a single component, or `nothing` if the component
has no state. Components extend this method for their own `Input` part.
"""
initial_state(part::InputPart, input::AbstractInput) = nothing

"""
    initialize!(part, input, state)

Finish initializing the `state`, once all its parts are created. Components
that need the complete state to initialize, extend this method for their own
`Input` part.
"""
initialize!(part::InputPart, input::AbstractInput, state::AbstractState) = nothing

"""
    initial_state(input)

The initial state for a model. The parts of the state are collected from all
components, after which `initialize!` is called for each component, in
dependency order.
"""
function initial_state(input::AbstractInput)
    state = map_components(p -> initial_state(p, input), input)
    foreach(p -> initialize!(p, input, state), parts(input))
    return state
end

end
# ~/~ end
