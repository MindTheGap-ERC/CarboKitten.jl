using CSV, DataFrames, GLMakie

# Step 1: Load the CSV data
# Replace "file.csv" with the path to your file
data = CSV.read("data/output/alcap2_sc.csv", DataFrame)

# Assuming column names: `time`, `thickness1`, `thickness2`, and `thickness3`
time = data[:,1]


litho1 = data[:,5]
litho2 = data[:,6]
litho3 = data[:,7]

sum_facies = litho1 .+ litho2 .+ litho3
sum_facies = vcat(sum,0)


fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1], xlabel = "Lithology", ylabel = "Time (Ma)")

# Define x-positions for lithologies (stacked columns)
x_left = 0.2  # Left boundary of each lithology column
x_right = 0.5  # Right boundary of each lithology column

# Plot each lithology as a filled rectangle (polygon)
for i in 1:length(time)-1
    boundary1 = zeros(Float64,length(time))
    boundary2 = zeros(Float64,length(time))
    boundary3 = zeros(Float64,length(time))

    boundary1[i+1] = sum_facies[i] + boundary1[i] # bot of Lithology 1
    boundary2[i+1] = sum_facies[i] .+ boundary1[i+1]  # bot of Lithology 2
    boundary3[i+1] = sum_facies[i] .+ boundary2[i+1] # bot of Lithology 3
    hspan!(ax,boundary1[i], boundary2[i+1],color = :green)
    hspan!(ax, boundary2[i], boundary3[i+1], color = :red)
    hspan!(ax, boundary3[i], boundary3[i+1],color = :blue)
end



# Customize plot aesthetics
ax.ylabel = "Thickness (m)"
hidedecorations!(ax, grid=false)  # Hide grid lines if needed
hidexdecorations!(ax)             # Hide x-axis tick labels if you don't want them
hideydecorations!(ax)             # Hide y-axis tick labels if you don't want them

# Save the figure as a high-resolution PNG
save("stratigraphic_column.png", fig, px_per_unit=2)