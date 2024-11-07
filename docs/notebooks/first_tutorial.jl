### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ ac7ec9d8-a70e-4b0e-be7b-705037273165
# Remove this once CarboKitten 0.2 is released
using Pkg; Pkg.activate("../../workenv")

# ╔═╡ 9c050dc8-b09f-404a-a699-d8ce88aa1def
using Revise

# ╔═╡ bcea7127-3c21-4c35-af42-3d2c71464409
using CarboKitten

# ╔═╡ d72d7e42-8392-44a0-a8b3-d59475be8dc7
using CarboKitten.Components.Common         # FIXME: get rid of this import

# ╔═╡ 325e3d04-2ff2-4c27-91bf-265820ac6763
using CarboKitten.Model.ALCAP2: run as run_model, Example, ALCAP2 as ALCAP

# ╔═╡ 4fe0f485-2db0-4685-b5b9-e9ba6009e4a6
using GLMakie

# ╔═╡ 5432b762-50fa-4ac4-97e6-0477236cb94a
using CarboKitten.Visualization: summary_plot

# ╔═╡ 9babc2b0-9c26-11ef-3459-0d113ec3c402
md"""
# Install Julia and Pluto

Please install Julia from the following webpage: [https://julialang.org/downloads/](https://julialang.org/downloads/).

We will use [Pluto](https://plutojl.org/) to do our tutorial. This is a notebook interface (similar to Jupyter) that is easy to use and has a strong focus on reproducibility.

Start the Julia REPL (read-eval-print loop), either from the start menu, or in a terminal, by running `julia`. You should see a colorful welcome message and a prompt for input:
"""

# ╔═╡ 17722d8b-baca-4f16-981f-1501c734a95f
md"""
```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.1 (2024-10-16)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```
"""

# ╔═╡ 22ec7e16-b8c0-414d-9700-52bf379e1051
md"""
You can install Pluto by running `using Pluto` and then answering `y` to the prompted question.

```juliarepl
julia> using Pluto
 │ Package Pluto not found, but a package named Pluto is available from a
 │ registry. 
 │ Install package?
 │   (@v1.11) pkg> add Pluto 
 └ (y/n/o) [y]: 
```

After a while you should see the following message:

```
┌ Info: 
│   Welcome to Pluto v0.20.3 🎈
│   Start a notebook server using:
│ 
│ julia> Pluto.run()
│ 
│   Have a look at the FAQ:
│   https://github.com/fonsp/Pluto.jl/wiki
└ 
```

You're good to go and run Pluto now!
"""

# ╔═╡ 3b7fef8b-efb9-467d-b6db-f7cfa132be69
md"""
## Install CarboKitten

In a Pluto notebook, Julia packages are installed by using them.
"""

# ╔═╡ cf61bc3f-a20a-45b7-a885-22b70075fc42
md"""
All packages used and their versions are stored inside the notebooks. When you run the notebook on a different computer, the same packages will be used.
"""

# ╔═╡ 68fac1d8-f402-429e-90a4-25fcfa188c2e
md"# A first example"

# ╔═╡ 74f4674f-dbea-44ad-8d54-9861b35139cd
run_model(Model{ALCAP}, Example.INPUT, "../../data/output/tutorial.h5")

# ╔═╡ e118a117-9a00-4589-906d-c31d2057bcef
summary_plot("../../data/output/tutorial.h5")

# ╔═╡ 8f883ea5-d90e-41e7-9809-3f170183a640


# ╔═╡ 69f42696-654e-4a82-b2de-a8a30fed22fb
md"The planview of the model is a 'chessboard' with 50 (along shore) × 100 (perpendicular to the shore) grids. The grid size is 150× 150m. With time, the produced carbonate would be added vertically to the 'chessboard'"

# ╔═╡ 25f1664e-a20e-4700-97f3-22e42bd6e3e7
md"`box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m)`"

# ╔═╡ 5426d044-7f3c-4ae6-b9c3-a82982f56f3c
md"This example assumes sea-level to be a sine curve"

# ╔═╡ abb330a5-b1d4-4b3d-8493-436ed96428c7
md" 'sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD)', where PERIOD is period of the curve = 0.2Myr in this case"

# ╔═╡ 08d1bf20-83b3-4275-a090-edeb207034a1
# ╠═╡ disabled = true
#=╠═╡
using GLMakie
  ╠═╡ =#

# ╔═╡ f6776dbe-8dc5-4c01-a034-cc47d0c34e8a
# ╠═╡ disabled = true
#=╠═╡
begin 
const PERIOD = 0.2
const AMPLITUDE = 4.0
t = collect(0:0.002:1)
sea_level= AMPLITUDE .* sin.(2π .* t ./ PERIOD)
fig1 = lines(t, sea_level, color = :blue, linewidth = 2)
fig1
end
  ╠═╡ =#

# ╔═╡ 5c90b464-858a-4a5d-a2c7-fda775057ef6
md"The initial topography is a ramp"

# ╔═╡ b4b13d2f-9920-4327-a9c2-893661744085
# ╠═╡ disabled = true
#=╠═╡
begin
x = collect(0:150:15000)
bedrock_elevation= -x / 300.0
fig2 = lines(x, bedrock_elevation, color = :black, linewidth = 2)
fig2
end
  ╠═╡ =#

# ╔═╡ a6d2e583-8469-43af-a50a-6b3dc7fbfe40
md"The subsidence_rate=50.0m / Myr, and insolation is a constant =400.0W/m^2"

# ╔═╡ 7829da51-0ff4-43c3-9a06-1fe7484c3a1a
md"Three biological facies are considered in this case, including: 
- T(Tropical) factory (max growth rate = 500.0m/Myr)
- M(Mounds) factory (max growth rate = 400.0m/Myr)
- C(Cool water) factory (max growth rate = 100.0m/Myr)
and their growth rates depends on the water depth (light penetration)"

# ╔═╡ 119c0993-4028-4d08-9e83-6d1d3683f93d
md"Let's first see how the cross-sections look like"

# ╔═╡ adda3ee8-c461-40e6-bf4e-ca676942ce3b
# ╠═╡ disabled = true
#=╠═╡
begin
using GLMakie
using CarboKitten.Visualization

save("docs/alcaps2_exmaple.png", summary_plot("data/output/alcap2.h5"))
# need to change the relative path after the simple example is done
end
  ╠═╡ =#

# ╔═╡ 7ff4d3e4-670f-4cab-9220-d1f8471efb50
md"The model would store the results in a HDF5 file, and automatically helps you extract age-depth model from locations from land towards sea"

# ╔═╡ cf65176e-1c3f-4059-897f-de645d72cd29
md"In this case, we extract the 10th, 30th, 50th and 70th grid from land (the position where red dashed lines are, the next figure)"

# ╔═╡ 5987c45f-bf4d-40f5-99b9-0ab9b2f25ada
# ╠═╡ disabled = true
#=╠═╡
begin
fig3 = lines(x, bedrock_elevation, color = :black, linewidth = 2)
sea_surface = AMPLITUDE .* sin.(2π .* x ./ PERIOD)
x_vlines = collect(10:20:70) .* 150
for x_pos in x_vlines
    vlines!(x_pos, color = :red, linestyle = :dash, linewidth = 1.5)
end
lines!(x,sea_surface,color = :blue, linewidth = 1)
fig3
end
  ╠═╡ =#

# ╔═╡ b4d995fd-dfc4-47c0-b428-d0dfc4a606ea
md"The blue line is sea-level, black line is initial topography and the red-dashed lines indicate the positions we are extracting data"

# ╔═╡ 48975f1c-55fb-4f72-a858-2f4cf301ac51
md"The data (i.e., age-depth model) are stored in `data/output`, with the name of the run. Herein, 'adm' = 'age-depth model', 'sc' = 'stratigraphic columns'"

# ╔═╡ 2d6569d2-f191-49ce-893c-7a4beaed0c96
# ╠═╡ disabled = true
#=╠═╡
begin
using DataFrames, CSV 
alcap2_example_adm = CSV.read("data/output/alcap2_adm.csv", DataFrame)
fig5 = Figure()
ax = Axis(fig5[1, 1],title = "Age-Depth Model Plots", xlabel = "Age (Myr)", ylabel = "Depth (m)")
time = alcap2_example_adm[1:end,1]
location = ["10","30","50","70"]
for (i, ycol) in enumerate([2, 3, 4, 5])
    fig5 = lines!(ax, time, alcap2_example_adm[!, ycol], label = location[i])
end
fig5[1,2] = Legend(fig5,ax)
fig5
end

  ╠═╡ =#

# ╔═╡ 4802af1d-8e60-4dc9-85b6-3f057be65336
md"# Small Tasks"

# ╔═╡ 6065af94-5dd4-485d-9c11-7a1360b1458c
md"Let's try to modify the sea-level curve. For example, add 2nd and 3rd order of sine curves o top of the existing sea-level curve. How this differs from the example? How about changing the Magnitude and Amplitudes, or gradient of initial topography? You might need to do some modifications on the input files (Shown later)"

# ╔═╡ 1e7c85db-7161-49f3-a54d-e6eef8bf799a
md"After this, you can try to a real sea-level curve"

# ╔═╡ 41d1bf67-1eb9-4272-8ffb-3c0c0ea7926e
md"The input file (consists of four parts) looks like:"

# ╔═╡ b9aa7141-ede6-47b2-a593-a601eddf0d61
md"first part: dependency. It also defines where you would like to store your output in `PATH` and `TAG`"

# ╔═╡ 3cc4223e-7ae0-4c42-aca1-14d86cb32a5d
# ╠═╡ disabled = true
#=╠═╡
begin
using Unitful
using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Model: ALCAP2 as ALCAP
using CarboKitten.Export: data_export, CSV

const m = u"m"
const Myr = u"Myr"

const PATH = "data/output"
const TAG = "alcap2"
end

  ╠═╡ =#

# ╔═╡ d5063ba9-9ebc-43b1-b262-c82272e6318e
md"second part: Facies defination"

# ╔═╡ d8d2d274-b63f-4364-9c25-5349ad7edb2c
# ╠═╡ disabled = true
#=╠═╡
begin
const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10000u"m"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5000u"m"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=7000u"m")
]
end
  ╠═╡ =#

# ╔═╡ 3a671070-332e-4411-9fed-a3945c04c830
md"the third part is the input values"

# ╔═╡ 7532c74b-6550-46c5-8ddd-bbc095c91780
# ╠═╡ disabled = true
#=╠═╡
begin
	const PERIOD = 0.2Myr
    const AMPLITUDE = 4.0m
	
	const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1),
    ca_interval=1,
    bedrock_elevation=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    facies=FACIES)
end
  ╠═╡ =#

# ╔═╡ 9acb27a7-cd33-4209-88bc-66268a1e1ee4
md"the fourth part is the data export"

# ╔═╡ 22170bcc-402a-4045-bbf6-a209c65907a8
# ╠═╡ disabled = true
#=╠═╡
begin
	function main()
    H5Writer.run(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSV(tuple.(10:20:70, 25),
          :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
          :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
          :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
          :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end
end
  ╠═╡ =#

# ╔═╡ 368cc3f8-80c0-4665-af56-33dd427bfa0e
md"run the model"

# ╔═╡ ec47c9f0-6cc0-41aa-b26f-d866e59f63e7
#=╠═╡
main()
  ╠═╡ =#

# ╔═╡ 3f2a50ac-1297-4a03-a252-bbc96b4f1137
md"You can also try your own sea-level curve. Herein filepath means the place where you store your sea-level curve. Make sure your sea-level curve has two columns."

# ╔═╡ 86aedf7b-fcaa-4eda-b1df-12cc38975e15
# ╠═╡ disabled = true
#=╠═╡
begin 
using Interpolations
function sealevel_curve(t,filepath)
    data = DataFrame(CSV.File(filepath))
	time = data[:,1] #assumes your first column is time
    sea_level = data[:,2] 
    x = linear_interpolation(time, sea_level)
    return x(t)
end 
end
  ╠═╡ =#

# ╔═╡ 123fa963-9a89-469e-a66f-b788d5575ff5
md"You can now replace the sea-level curve with: `sea_level = t -> sealevel_curve(t,filepath)` "

# ╔═╡ Cell order:
# ╠═ac7ec9d8-a70e-4b0e-be7b-705037273165
# ╠═9c050dc8-b09f-404a-a699-d8ce88aa1def
# ╟─9babc2b0-9c26-11ef-3459-0d113ec3c402
# ╟─17722d8b-baca-4f16-981f-1501c734a95f
# ╟─22ec7e16-b8c0-414d-9700-52bf379e1051
# ╟─3b7fef8b-efb9-467d-b6db-f7cfa132be69
# ╠═bcea7127-3c21-4c35-af42-3d2c71464409
# ╟─cf61bc3f-a20a-45b7-a885-22b70075fc42
# ╠═68fac1d8-f402-429e-90a4-25fcfa188c2e
# ╠═d72d7e42-8392-44a0-a8b3-d59475be8dc7
# ╠═325e3d04-2ff2-4c27-91bf-265820ac6763
# ╠═74f4674f-dbea-44ad-8d54-9861b35139cd
# ╠═4fe0f485-2db0-4685-b5b9-e9ba6009e4a6
# ╠═5432b762-50fa-4ac4-97e6-0477236cb94a
# ╠═e118a117-9a00-4589-906d-c31d2057bcef
# ╠═8f883ea5-d90e-41e7-9809-3f170183a640
# ╠═69f42696-654e-4a82-b2de-a8a30fed22fb
# ╠═25f1664e-a20e-4700-97f3-22e42bd6e3e7
# ╠═5426d044-7f3c-4ae6-b9c3-a82982f56f3c
# ╠═abb330a5-b1d4-4b3d-8493-436ed96428c7
# ╠═08d1bf20-83b3-4275-a090-edeb207034a1
# ╠═f6776dbe-8dc5-4c01-a034-cc47d0c34e8a
# ╠═5c90b464-858a-4a5d-a2c7-fda775057ef6
# ╠═b4b13d2f-9920-4327-a9c2-893661744085
# ╠═a6d2e583-8469-43af-a50a-6b3dc7fbfe40
# ╠═7829da51-0ff4-43c3-9a06-1fe7484c3a1a
# ╠═119c0993-4028-4d08-9e83-6d1d3683f93d
# ╠═adda3ee8-c461-40e6-bf4e-ca676942ce3b
# ╠═7ff4d3e4-670f-4cab-9220-d1f8471efb50
# ╠═cf65176e-1c3f-4059-897f-de645d72cd29
# ╠═5987c45f-bf4d-40f5-99b9-0ab9b2f25ada
# ╠═b4d995fd-dfc4-47c0-b428-d0dfc4a606ea
# ╠═48975f1c-55fb-4f72-a858-2f4cf301ac51
# ╠═2d6569d2-f191-49ce-893c-7a4beaed0c96
# ╠═4802af1d-8e60-4dc9-85b6-3f057be65336
# ╠═6065af94-5dd4-485d-9c11-7a1360b1458c
# ╠═1e7c85db-7161-49f3-a54d-e6eef8bf799a
# ╠═41d1bf67-1eb9-4272-8ffb-3c0c0ea7926e
# ╠═b9aa7141-ede6-47b2-a593-a601eddf0d61
# ╠═3cc4223e-7ae0-4c42-aca1-14d86cb32a5d
# ╠═d5063ba9-9ebc-43b1-b262-c82272e6318e
# ╠═d8d2d274-b63f-4364-9c25-5349ad7edb2c
# ╠═3a671070-332e-4411-9fed-a3945c04c830
# ╠═7532c74b-6550-46c5-8ddd-bbc095c91780
# ╠═9acb27a7-cd33-4209-88bc-66268a1e1ee4
# ╠═22170bcc-402a-4045-bbf6-a209c65907a8
# ╠═368cc3f8-80c0-4665-af56-33dd427bfa0e
# ╠═ec47c9f0-6cc0-41aa-b26f-d866e59f63e7
# ╠═3f2a50ac-1297-4a03-a252-bbc96b4f1137
# ╠═86aedf7b-fcaa-4eda-b1df-12cc38975e15
# ╠═123fa963-9a89-469e-a66f-b788d5575ff5
