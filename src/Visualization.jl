# ~/~ begin <<docs/src/visualization.md#src/Visualization.jl>>[init]
module Visualization
export sediment_profile!, sediment_profile, wheeler_diagram!, wheeler_diagram, production_curve!, production_curve

print_instructions() = print("This is an extension and only becomes available when you import {Cairo,GL,WGL}Makie before using this.\n")

sediment_profile!(args...) = print_instructions()
sediment_profile(args...) = print_instructions()
wheeler_diagram!(args...) = print_instructions()
wheeler_diagram(args...) = print_instructions()
production_curve(args...) = print_instructions()
production_curve!(args...) = print_instructions()
end  # module
# ~/~ end
