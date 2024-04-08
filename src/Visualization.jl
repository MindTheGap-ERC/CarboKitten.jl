# ~/~ begin <<docs/src/ca-with-production.md#src/Visualization.jl>>[init]
module Visualization
    export plot_crosssection, plot_facies_production

    print_instructions() = print("You'll need to import both Makie and GeometryBasics before using this.\n")

    function plot_facies_production(args...)
        print_instructions()
    end

    function plot_crosssection(args...)
        print_instructions()
    end
end  # module
# ~/~ end