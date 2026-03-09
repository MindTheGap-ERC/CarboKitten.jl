# ~/~ begin <<docs/src/components/cellular-automata.md#examples/model/alcap/feedback-plot.jl>>[init]
module PlotFeedback
    using GLMakie
    using CarboKitten.Visualization
    using CarboKitten.Export: read_volume, read_slice

    function main()
        GLMakie.activate!()
        fig = Figure(size=(1000, 800))

        header, topography = read_volume("data/output/ca-wo-feedback.h5", :topography)
        _, profile = read_slice("data/output/ca-wo-feedback.h5", :profile)

        ax1 = Axis(fig[1, 1:2])
        sediment_profile!(ax1, header, profile)
        ax1.title = "without feedback"
        
        ax1_wh1 = Axis(fig[3, 1])
        ax1_wh2 = Axis(fig[3, 2])
        sa, ft = wheeler_diagram!(ax1_wh1, ax1_wh2, header, profile)
        Colorbar(fig[2, 1], sa; vertical=false, label="sediment accumulation [m/Myr]")
        Colorbar(fig[2, 2], ft; vertical=false, ticks=1:3, label="dominant facies")

        header, topography = read_volume("data/output/ca-feedback.h5", :topography)
        _, profile = read_slice("data/output/ca-feedback.h5", :profile)

        ax2 = Axis(fig[1, 3:4])
        sediment_profile!(ax2, header, profile)
        ax2.title = "with feedback"
        
        ax2_wh1 = Axis(fig[3, 3])
        ax2_wh2 = Axis(fig[3, 4])
        sa, ft = wheeler_diagram!(ax2_wh1, ax2_wh2, header, profile)

        Colorbar(fig[2, 3], sa; vertical=false, label="sediment accumulation [m/Myr]")
        Colorbar(fig[2, 4], ft; vertical=false, ticks=1:3, label="dominant facies")
        fig
    end
end

save("docs/src/_fig/ca-feedback.png", PlotFeedback.main())
# ~/~ end
