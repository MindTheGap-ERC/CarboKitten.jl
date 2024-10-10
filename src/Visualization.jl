# ~/~ begin <<docs/src/visualization.md#src/Visualization.jl>>[init]
module Visualization
export sediment_profile!, sediment_profile, wheeler_diagram!, wheeler_diagram, production_curve!, production_curve

function print_instructions(func_name, args)
    println("Called `$(func_name)` with args `$(typeof.(args))`")
    println("This is an extension and only becomes available when you import {Cairo,GL,WGL}Makie before using this.")
end

sediment_profile!(args...) = print_instructions("sediment_profile!", args)
sediment_profile(args...) = print_instructions("sediment_profile", args)
wheeler_diagram!(args...) = print_instructions("wheeler_diagram!", args)
wheeler_diagram(args...) = print_instructions("wheeler_diagram", args)
production_curve(args...) = print_instructions("production_curve", args)
production_curve!(args...) = print_instructions("production_curve!", args)
stratigraphic_column!(args...) = print_instructions("production_curve!", args)
age_depth_model!(args...) = print_instructions("age_depth_model!", args)

end  # module
# ~/~ end