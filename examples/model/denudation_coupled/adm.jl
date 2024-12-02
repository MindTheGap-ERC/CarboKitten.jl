using CSV, DataFrames, GLMakie

dis_adm = CSV.read("data/output/dissolution_adm.csv", DataFrame)
phys_adm = CSV.read("data/output/physical_adm.csv", DataFrame)
emp_adm = CSV.read("data/output/empirical_adm.csv", DataFrame)

time = dis_adm[1:end,1]
dis_adm_shallow = dis_adm[1:end,2]
phys_adm_shallow = phys_adm[1:end,2]
emp_adm_shallow = emp_adm[1:end,2]
adm_total = hcat(time,dis_adm_shallow,phys_adm_shallow,emp_adm_shallow)
fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1],title = "Age-Depth Model Plots", xlabel = "Age (Myr)", ylabel = "Depth (m)")
mode = ["Dissolution","Physical Erosion", "Empirical Denudation"]

for i in 1:3
    lines!(ax, time, adm_total[:,i+1], label = mode[i])
end

fig