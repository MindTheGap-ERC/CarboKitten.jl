# ~/~ begin <<docs/src/visualization.md#src/Visualization.jl>>[init]
module Visualization
    export plot_crosssection, plot_facies_production, sediment_profile

    print_instructions() = print("This is an extension and only becomes available when you import {Cairo,GL,WGL}Makie before using this.\n")

    plot_facies_production(args...) = print_instructions()
    sediment_profile!(args...) = print_instructions()
    sediment_profile(args...) = print_instructions()
end  # module
# ~/~ end
