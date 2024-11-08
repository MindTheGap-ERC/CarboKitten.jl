using DataFrames, CSV, GLMakie
alcap2_example_adm = CSV.read("data/output/alcap2_adm.csv", DataFrame)
fig = Figure()
ax = Axis(fig[1, 1],title = "Age-Depth Model Plots", xlabel = "Age (Myr)", ylabel = "Depth (m)")
time = alcap2_example_adm[1:end,1]
location = ["10","30","50","70"]
for (i, ycol) in enumerate([2, 3, 4, 5])
    fig = lines!(ax, time, alcap2_example_adm[!, ycol], label = location[i])
end
fig[1,2] = Legend(fig,ax)
fig
save("docs/src/_fig/adm.png", fig) 