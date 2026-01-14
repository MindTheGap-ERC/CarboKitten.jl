# ~/~ begin <<docs/src/boxes.md#src/Config.jl>>[init]
module Config

export TimeProperties

using ..BoundaryTrait
using ..Vectors

using ..Boxes: Box, axes
export Box, axes

using Unitful
using Unitful.DefaultSymbols
using ..CarboKitten: TimeProperties

# ~/~ begin <<docs/src/boxes.md#config-types>>[init]

# ~/~ end

end
# ~/~ end
