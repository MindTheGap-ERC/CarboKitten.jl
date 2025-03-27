### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 66816e88-a59c-4ce9-93fa-2a9959b740f1
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 6941c10a-097f-11f0-1abd-2bd1f922bf6e
module Script

using Unitful
using CarboKitten
using CarboKitten.Export: data_export, CSV
using CarboKitten.Transport.Solvers: forward_euler, runge_kutta_4

v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

wave_velocity(v_max) = w -> let (v, s) = v_prof(v_max, 10.0u"m", w)
		((v, 0.0u"m/Myr"), (s, 0.0u"1/Myr"))
	end

const PATH = "../data/output"
const TAG = "ot-example"

const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr",
		wave_velocity=wave_velocity(-0.01u"m/yr")),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10.0u"m/yr",
		wave_velocity=wave_velocity(-1.0u"m/yr")),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5.0u"m/yr",
		wave_velocity=wave_velocity(-0.2u"m/yr"))
]

const PERIOD = 0.2u"Myr"
const AMPLITUDE = 4.0u"m"
const BOX = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")
const INPUT = ALCAP.Input(
    tag="$TAG",
    box=BOX,
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000,
        write_interval=1),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    transport_solver=runge_kutta_4(typeof(1.0u"m"), BOX),
	transport_substeps=10,
    facies=FACIES)

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSV(tuple.(10:20:70, 25),
            :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
            :water_depth => "$(PATH)/$(TAG)_wd.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end


# ╔═╡ 8e1bf0aa-82c9-4e6e-98c2-fdf05e4ce816
using GLMakie

# ╔═╡ 42ae3cbb-b02a-4f88-926d-67e337299412
using CarboKitten.Visualization: summary_plot

# ╔═╡ f755055d-ee7f-4951-ac13-34f8a13f9396
module Island

using Unitful
using CarboKitten
using CarboKitten.Export: data_export, CSV
using CarboKitten.Transport.Solvers: forward_euler, runge_kutta_4

v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

wave_velocity(v_max) = w -> let (v, s) = v_prof(v_max, 10.0u"m", w)
		((v, 0.0u"m/Myr"), (s, 0.0u"1/Myr"))
	end

initial_topography(x, y) = -sqrt((x - 7.5u"km")^2 + (y - 7.5u"km")^2) / 150.0

const PATH = "../data/output"
const TAG = "ot-island"

const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10.0u"m/yr",
		wave_velocity=wave_velocity(-0.01u"m/yr")),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5.0u"m/yr",
		wave_velocity=wave_velocity(-1.0u"m/yr")),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5.0u"m/yr",
		wave_velocity=wave_velocity(-1.0u"m/yr"))
]

const PERIOD = 0.2u"Myr"
const AMPLITUDE = 4.0u"m"
const BOX = Box{Periodic{2}}(grid_size=(50, 50), phys_scale=300.0u"m")
const INPUT = ALCAP.Input(
    tag="$TAG",
    box=BOX,
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000,
        write_interval=1),
    ca_interval=1,
    initial_topography=initial_topography,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    transport_solver=forward_euler, # runge_kutta_4(typeof(1.0u"m"), BOX),
	transport_substeps=5,
    facies=FACIES)

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")
end

end


# ╔═╡ 5fb3fffe-3ee3-47c7-abb8-efd166a0013a
using HDF5

# ╔═╡ a02c8a4b-e265-4d73-8f1a-8ca988d8c09e
using CarboKitten.Export: read_header

# ╔═╡ 4c88a358-c083-4882-8ce4-b6ccc1c6f167
using Unitful

# ╔═╡ 4331ed00-8489-4a34-9e19-984798ac815c
# Script.main()

# ╔═╡ e857236c-a7ce-4f74-9280-d4dc03ccaba8
summary_plot("../data/output/ot-example.h5")

# ╔═╡ 71c098a6-864b-45b5-b365-6d3eed5baec2
Island.main()

# ╔═╡ 3ec2198e-fbae-46ce-92f7-eb0ce3830885
fig = summary_plot("../data/output/ot-island.h5")

# ╔═╡ 56806154-db7d-4ae6-beff-bc6f963d8e64
save("ot-island.png", fig)

# ╔═╡ 0e770c25-a0a8-4edc-89bc-8728c57ed5d9
h5open("../data/output/ot-island.h5") do fid
	h = read_header(fid)
	s = fid["sediment_height"][:, :, end] * u"m"
	t = h.initial_topography .+ s .- (h.axes.t[end] * h.subsidence_rate)

	fig = Figure()
	ax = Axis(fig[1, 1])
	hm = heatmap!(ax, h.axes.x, h.axes.y, t / u"m", colorrange=(-2.0, 1.0), colormap=:curl)
	Colorbar(fig[1, 2], hm)

	fig
end

# ╔═╡ Cell order:
# ╠═66816e88-a59c-4ce9-93fa-2a9959b740f1
# ╠═6941c10a-097f-11f0-1abd-2bd1f922bf6e
# ╠═4331ed00-8489-4a34-9e19-984798ac815c
# ╠═8e1bf0aa-82c9-4e6e-98c2-fdf05e4ce816
# ╠═42ae3cbb-b02a-4f88-926d-67e337299412
# ╠═e857236c-a7ce-4f74-9280-d4dc03ccaba8
# ╠═f755055d-ee7f-4951-ac13-34f8a13f9396
# ╠═71c098a6-864b-45b5-b365-6d3eed5baec2
# ╠═3ec2198e-fbae-46ce-92f7-eb0ce3830885
# ╠═56806154-db7d-4ae6-beff-bc6f963d8e64
# ╠═5fb3fffe-3ee3-47c7-abb8-efd166a0013a
# ╠═a02c8a4b-e265-4d73-8f1a-8ca988d8c09e
# ╠═4c88a358-c083-4882-8ce4-b6ccc1c6f167
# ╠═0e770c25-a0a8-4edc-89bc-8728c57ed5d9
