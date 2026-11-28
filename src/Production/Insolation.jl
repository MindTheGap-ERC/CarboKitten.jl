# ~/~ begin <<docs/src/production/production.md#src/Production/Insolation.jl>>[init]
module Insolation

using Unitful
using Interpolations: linear_interpolation, Flat

import ..Abstract: AbstractInput
import ...CarboKitten: time_axis

"""
    insolation_curve(input) -> function(time) -> insolation

Build a closure mapping time to insolation from the input specification.
Handles constant (`Quantity`), tabular (`AbstractVector`), and functional
insolation inputs.
"""
function insolation_curve(input::AbstractInput)
    insolation_param = input.insolation
    if insolation_param isa Quantity
        return _ -> insolation_param
    elseif insolation_param isa AbstractVector
        t_axis = time_axis(input.time)
        length(insolation_param) == length(t_axis) || throw(ArgumentError(
            "tabular insolation needs one value per time axis point " *
            "($(length(t_axis)) for $(input.time.steps) steps), got $(length(insolation_param))"))
        t_vals = ustrip.(u"Myr", t_axis)
        I_vals = ustrip.(u"W/m^2", insolation_param)
        itp = linear_interpolation(t_vals, I_vals, extrapolation_bc=Flat())
        return t -> itp(ustrip(u"Myr", t)) * u"W/m^2"
    else
        return t -> insolation_param(t)
    end
end

end
# ~/~ end
