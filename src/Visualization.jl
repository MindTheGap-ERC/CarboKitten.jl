# ~/~ begin <<docs/src/visualization.md#src/Visualization.jl>>[init]
module Visualization
export sediment_profile!, sediment_profile,
       wheeler_diagram!, wheeler_diagram,
       production_curve!, production_curve,
       glamour_view!, glamour_view,   # ← ADD THIS
       coeval_lines!, profile_plot!, summary_plot,
       dominant_facies!, sediment_accumulation!,
       summary_plot_from_state, profile_plot_coarse!, extract_column, map_view, wheeler_column, fence_plot, sediment_profile_generic!, build_vertical_coordinates, resample_to_regular_grid, stratigraphic_column_layers!, facies_colormap



function print_instructions(func_name, args)
    println("Called `$(func_name)` with args `$(typeof.(args))`")
    println("This is an extension and only becomes available when you import {Cairo,GL,WGL}Makie before using this.")
end

function build_vertical_coordinates end

function resample_to_regular_grid end


function profile_plot! end

# profile_plot!(args...; kwargs...) = print_instructions("profile_plot!", args)
sediment_accumulation!(args...) = print_instructions("sediment_accumulation!", args)
dominant_facies!(args...) = print_instructions("dominant_facies!", args)


coeval_lines!(args...) = print_instructions("coeval_lines!", args)


sediment_profile!(args...) = print_instructions("sediment_profile!", args)
sediment_profile(args...) = print_instructions("sediment_profile", args)


wheeler_diagram!(args...) = print_instructions("wheeler_diagram!", args)


wheeler_diagram(args...) = print_instructions("wheeler_diagram", args)
production_curve(args...) = print_instructions("production_curve", args)


production_curve!(args...) = print_instructions("production_curve!", args)


stratigraphic_column!(args...) = print_instructions("production_curve!", args)


glamour_view!(args...) = print_instructions("glamour_view!", args)


glamour_view(args...) = print_instructions("glamour_view", args)


summary_plot(args...) = print_instructions("summary_plot", args)


function summary_plot_from_state end
function profile_plot_coarse! end


function fence_plot end


function extract_column end


function map_view end


function wheeler_column end


function sediment_profile_generic! end

function stratigraphic_column_layers! end

function facies_colormap end

end  # module
# ~/~ end
