# ~/~ begin <<docs/src/components/boxes.md#src/Config.jl>>[init]
module Config

import ..CarboKitten

Base.@deprecate_binding TimeProperties CarboKitten.TimeProperties
Base.@deprecate_binding Box CarboKitten.Box
Base.@deprecate_binding axes CarboKitten.box_axes

end
# ~/~ end
