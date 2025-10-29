# ~/~ begin <<docs/src/getting-started.md#examples/getting-started/plot_alcap.jl>>[init]
using GLMakie
using CarboKitten.Visualization

# Activate the GLMakie backend
GLMakie.activate!()

# Generate a summary plot
fig = summary_plot("data/output/first-run.h5")

# Display the plot
display(fig)

# Optionally save to file
save("my-first-model.png", fig)
# ~/~ end
