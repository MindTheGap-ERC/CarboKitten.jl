using CSV, DataFrames, GLMakie

data = CSV.read("data/output/alcap2_sc.csv", DataFrame)
time = data[:,1]
litho1 = data[:,5]
litho2 = data[:,6]
litho3 = data[:,7]

sum_facies = litho1 .+ litho2 .+ litho3
sum_facies = vcat(sum_facies,0)


fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1],  ylabel = "Thickness (m)")


x_left = 0.2  
x_right = 0.25  
boundary1 = zeros(Float64,length(time))
boundary2 = zeros(Float64,length(time))
boundary3 = zeros(Float64,length(time))


for i in 1:length(time)-1

    boundary1[i+1] = sum_facies[i] + boundary1[i] # bot of Lithology 1
    boundary2[i+1] = sum_facies[i] .+ boundary1[i+1]  # bot of Lithology 2
    boundary3[i+1] = sum_facies[i] .+ boundary2[i+1] # bot of Lithology 3
    hspan!(ax,boundary1[i], boundary1[i+1],color = Makie.wong_colors()[1])
    hspan!(ax, boundary2[i], boundary2[i+1], color = Makie.wong_colors()[2])
    hspan!(ax, boundary3[i], boundary3[i+1],color = Makie.wong_colors()[3])
end

save("docs/src/_fig/stratigraphic_column.png", fig, px_per_unit=2)