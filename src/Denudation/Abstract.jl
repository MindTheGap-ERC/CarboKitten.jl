# ~/~ begin <<docs/src/denudation/denudation.md#src/Denudation/Abstract.jl>>[init]
module Abstract

using ...BoundaryTrait: Boundary
using ...Boxes: Box
using ...SedimentStack: peek_sediment

using Unitful

abstract type DenudationType end

"""
    denudation(box, param, state)


Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.

Returns denudation mass in units of meters.
"""
function denudation(input)

    function (state, water_depth, slope)
        if denudation(input.box, input.denudation, water_depth, slope, input.facies, state) !== nothing
        return denudation(input.box, input.denudation, water_depth, slope, input.facies, state) .* input.time.Δt
        else
        return nothing
        end
    end
end

"""
    denudation(box::Box, param::DenudationType, water_depth, slope, facies)

Computes the amount of denudation. This function is called on a pixel by pixel basis, so all arguments can be assumed to be scalar. The `param` argument should be of a subtype of `DenudationType` containing all the input parameters for this specific denudation model.
"""
function denudation(box::Box, param::DenudationType, water_depth, slope, facies, state)
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

Takes `state`, `water_depth` in meters and `denudation_mass` as a 3D array (facies, x and y coordinates) in units of meters.
"""
function redistribution(input)
    function (state, water_depth, denudation_mass)
        return redistribution(input.box, input.denudation, denudation_mass, water_depth)
    end
end

function redistribution(box::Box, param::DenudationType, denudation_mass, water_depth)
    error("Abstract `redistribution` function called.")
end

# not sure this is the right place for this, but it's common to multiple modes of denudation
function dominant_facies(state, i::CartesianIndex, peek_depth::Float64)
    # look at top of the sediment buffer column, first two cells of buffer
    buffer_facies = peek_sediment(state.sediment_buffer[:,:,i[1],i[2]], peek_depth)
    max_f = findmax(buffer_facies)

    # we shouldn't be calling this function with an empty sediment buffer
    if max_f[1]==0.0 || isnan(max_f[1])
        @error "maximum facies value is $(max_f[1]), cannot find dominant facies if there's no sediment in buffer"
    else 
        return max_f[2]
    end
end

end  # module
# ~/~ end
