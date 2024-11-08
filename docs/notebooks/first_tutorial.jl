### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ ac7ec9d8-a70e-4b0e-be7b-705037273165
# Remove this once CarboKitten 0.2 is released
using Pkg; Pkg.activate("../../workenv")

# ╔═╡ 9c050dc8-b09f-404a-a699-d8ce88aa1def
using Revise

# ╔═╡ ebcb24b3-1ea3-49a8-a6d8-bf1f2cee657e
using PlutoUI

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

# ╔═╡ 51f2b3db-7b2d-4d62-8c3a-6608c01bffb7
using Unitful

# ╔═╡ b6ae4907-ac04-4a37-975f-42d85d572f91
using CarboKitten.Boxes: Box, Coast

# ╔═╡ 37cacf0c-c33a-4a11-b47e-fe2dc24d4eff
using CarboKitten.Config: TimeProperties

# ╔═╡ 90fdaca0-2efa-4230-8af3-177a6d94639c
using CarboKitten.Components.TimeIntegration: write_times

# ╔═╡ d0bfae1a-deef-4625-a242-0a9b899bf83d
using CarboKitten.Model.ALCAP2: Facies

# ╔═╡ 3723832f-344c-4471-acb0-cef7d4e5ca94
using CarboKitten.Export: data_export, CSV

# ╔═╡ 31e7c759-f980-4618-be90-892865751e58
using DataFrames

# ╔═╡ 61ae751b-0c05-4e96-8be9-3f85cb6afc51
using CSV: read as read_csv

# ╔═╡ 3c4cef70-df77-46ba-b623-fd46b5500e51
TableOfContents()

# ╔═╡ 0ce8de55-3304-431d-a2aa-110b46a25c9b
md"""
This CarboKitten tutorial is aimed at people that are completely new to Julia.
"""

# ╔═╡ 9babc2b0-9c26-11ef-3459-0d113ec3c402
md"""
# Installing

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
md"## A first example"

# ╔═╡ 9aafba01-fc4c-4dc1-85b6-9f33a4cfc77a
md"""Please make sure to set the output directory to a convenient place. If you downloaded this notebook to an empty directory, using `"."` would be a good choice.
"""

# ╔═╡ b3b271cb-143f-44ba-a593-80b9e6c96392
OUTPUTDIR = "../../data/output"

# ╔═╡ 1d5cf6cc-745d-4a5a-80ae-b1b6c057af0b
md"""
!!! tip "HDF5"
	Output is written to HDF5 files, a broadly used and supported data format for binary scientific data.

	After the simulation is done, we can read the HDF5 file for further analysis and visualization.
"""

# ╔═╡ c974cb9a-e1c0-4402-880c-7990d217da89
md"""
!!! info "Example: String interpolation"
	Notice the `"$(...)"` syntax: this is called string interpolation: the contents of the `OUTPUTDIR` variable are formatted into the resulting string. Try it out with a few examples (inside a `let` block so the variables don't leak):

	```julia
	let
		name = "Peter"
		"Hello, $(name)!"
	end
	```

	What happens if you change `name` into a number?
"""

# ╔═╡ 316b049d-698d-4d5e-9c18-73701ef8b492
md"""
!!! tip "Disabled cells"
    The cell below is disabled because it will take a minute or so to run. Click the cell-menu at the top right of the cell to enable it.
"""

# ╔═╡ 74f4674f-dbea-44ad-8d54-9861b35139cd
# ╠═╡ disabled = true
#=╠═╡
run_model(Model{ALCAP}, Example.INPUT, "$(OUTPUTDIR)/example.h5")
  ╠═╡ =#

# ╔═╡ 66e30189-ae72-4ec1-b7bd-1136ddfce2ee
md"""
!!! tip "Makie"
	We will be using [Makie](https://makie.org) to do our plotting.
"""

# ╔═╡ e118a117-9a00-4589-906d-c31d2057bcef
# ╠═╡ disabled = true
#=╠═╡
summary_plot("$(OUTPUTDIR)/example.h5")
  ╠═╡ =#

# ╔═╡ 8f883ea5-d90e-41e7-9809-3f170183a640
md"""
# A second example

In our second example we walk through the entire configuration. We start by defining the shape of our box.

## Simulation box
"""

# ╔═╡ 342490a7-685f-47f4-9889-a56927713f85
box = Box{Coast}(grid_size = (50, 50), phys_scale = 300.0u"m")

# ╔═╡ 5e914211-fd1f-45e0-9ecd-442d819684a2
md"""
There is a lot to unpack here: `Coast` indicates a box topology that is periodic in the y-direction and mirrored in the x-direction, which is ideal for simulating a strip of coast-line. Other choices could be `Periodic{2}` or `Reflected{2}`. The `grid_size` argument is the size of our box in pixels, and `phys_scale` is the size of one pixel with dimensions of length.
"""

# ╔═╡ a39a9476-6a7f-435a-b882-36ecf618369e
md"""
!!! info "Exercise: Units"
	To indicate units, we use the [Unitful package](https://painterqubits.github.io/Unitful.jl/stable/). This allows us to be free in our choice of units, as long as dimensions are in order.

	Try adding a quantity in light years to one in kilometers:

	```julia
	1.0u"yr" * 3e8u"m/s" + 150.0e6u"km"
	```

	Try adding a quantity in meters to one in years:

	```julia
	2.0u"m" + 3.0u"yr"
	```
"""

# ╔═╡ 620627e8-4bf6-400d-ab6f-9a7a9a6ec040
md"""
## Simulation time
Next, we define our simulation time and time step."""

# ╔═╡ 281df21c-aa0a-41fb-a20d-45f7b4eb3532
time = TimeProperties(
	Δt = 500u"yr",
	steps = 2000
)

# ╔═╡ a2689006-cb86-4956-999a-7adea546abdd
md"""Notice that the Δt property was automatically coneverted to Myr, which is the unit that we use internally.
"""

# ╔═╡ a81dba9b-8f12-4b72-9ea8-993d6c5e501b
md"""
## Sea level

Sea level is given as a function. The function should take quantities of time and return a quantity of length.

!!! tip "Functions in Julia"
	There are several ways to define functions in Julia. The most readable is as follows:

	```julia
	function some_function_name(arg1, arg2, etc...)
		# computations
		...
		return result
	end
	```

Here, we define the `sea_level` as a sine wave with fixed amplitude and period.
"""

# ╔═╡ 7d72a3bb-5ed8-4ebf-a45d-59e1bb267e9a
function sea_level(t)
	amplitude = 4.0u"m"
	period = 0.2u"Myr"
	
	return amplitude * sin(2π * t / period)
end

# ╔═╡ 30245436-9731-4cae-a150-fcba4360527e
let
	t = write_times(time)
	sl = sea_level.(t)
	lines(t, sl)
end

# ╔═╡ f6776dbe-8dc5-4c01-a034-cc47d0c34e8a
md"""
!!! info "Exercise: expand the `sea_level` function"
	Try to add a second mode of oscillation into the sea level curve. Choose your own amplitude and period.
"""

# ╔═╡ 5c90b464-858a-4a5d-a2c7-fda775057ef6
md"""
## Initial topography
In CarboKitten, the initial topography is also known as *bedrock elevation*. This quantity needs to be given as a function of both `x` and `y`. Here our initial topography is a ramp.
"""

# ╔═╡ 68f2dd79-6777-4267-ae67-5167f764b7b9
function bedrock_elevation(x, y)
	return -x / 300.0
end

# ╔═╡ 2e168646-5b9e-437b-b0b2-d637a4beb577
md"""
## Facies definitions

Three biological facies are considered in this case, including:

- T(Tropical) factory (max growth rate = 500.0m/Myr)
- M(Mounds) factory (max growth rate = 400.0m/Myr)
- C(Cool water) factory (max growth rate = 100.0m/Myr)

Their growth rates depend on the water depth, as parametrized by the Bosscher & Schlager 1992 model,

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right),$$

where $g_m$ is the maximum growth rate, $I_0$ is the insolation, $I_k$ the saturation intensity, $k$ the extinction coefficient, and $w$ the water depth.

Because our model also includes diffusive transport, we also specify a diffusion coefficient for each facies type. The units are a bit weird, but the intuition should be: larger values for smaller particle size.

When we inpsect the displayed `Facies` values, we see that also the default viability and activation ranges are set to govern the cellular automaton, following Burgess 2013.
"""

# ╔═╡ 68ab2b1c-67fb-48f1-ac81-c7efac04d96a
facies = [
	Facies(
		maximum_growth_rate = 500.0u"m/Myr",
		extinction_coefficient = 0.8u"m^-1",
		saturation_intensity = 60.0u"W/m^2",
		diffusion_coefficient = 500u"m"
	),
	Facies(
		maximum_growth_rate = 400.0u"m/Myr",
		extinction_coefficient = 0.1u"m^-1",
		saturation_intensity = 60.0u"W/m^2",
		diffusion_coefficient = 5000u"m"
	),
	Facies(
		maximum_growth_rate = 100.0u"m/Myr",
		extinction_coefficient = 0.005u"m^-1",
		saturation_intensity = 60.0u"W/m^2",
		diffusion_coefficient = 3000u"m"
	)]

# ╔═╡ 43405e49-2739-4e79-8204-e489db6c1fd5
md"""
## Other parameters

The subsidence rate is set to a constant 50m/Myr and insolation to 400W/m^2. The other parameters are related to the transport model.
"""

# ╔═╡ ae3a74ea-f6b5-442c-973b-1cec48627968
input = ALCAP.Input(
	tag = "tutorial",
	
	time = time,
	box = box,
	facies = facies,
	sea_level = sea_level,
	bedrock_elevation = bedrock_elevation,
	
	subsidence_rate = 50.0u"m/Myr",
	insolation = 400.0u"W/m^2",
	
	sediment_buffer_size = 50,
	depositional_resolution = 0.5u"m",
	
	disintegration_rate = 50.0u"m/Myr"
)

# ╔═╡ 8946d268-4407-4fe4-86ae-67b3a37b34be
# ╠═╡ disabled = true
#=╠═╡
run_model(Model{ALCAP}, input, "$(OUTPUTDIR)/tutorial2.h5")
  ╠═╡ =#

# ╔═╡ 0fd7ea33-ff06-4355-9278-125c8ed66df4
# ╠═╡ disabled = true
#=╠═╡
summary_plot("$(OUTPUTDIR)/tutorial2.h5")
  ╠═╡ =#

# ╔═╡ a3e5f420-a59d-4725-8f8f-e5b8f06987db
md"""
# Extracting CSV data

HDF5 is a very versatile data format, and is readable from every major computational platform (e.g. MatLab, Python, R). However, sometimes you may want to process your output data further using existing tools that read CSV data.
"""

# ╔═╡ ec4f0021-12af-4f8c-8acb-970e6820d2a4
export_locations = [(10, 25), (25, 25), (40, 25)]

# ╔═╡ 7f05b817-f933-4280-b2ed-ae318a535123
data_export(CSV(export_locations,
	:sediment_accumulation_curve => "$(OUTPUTDIR)/tutorial_sac.csv",
	:age_depth_model => "$(OUTPUTDIR)/tutorial_adm.csv"),
	"$(OUTPUTDIR)/tutorial.h5")

# ╔═╡ 2a24237e-5c1f-47e1-8f33-cca3ef563930
md"""
!!! info "Exercise: Meta data"
	Next to `:sediment_accumulation_curve` and `:age_depth_model`, we have implemented functions to extract `:stratigraphic_column` and `:metadata`.

	The `:metadata` target is special, since it doesn't write to CSV but to TOML. Add the `:metadata` target to the export list and inspect the result.
"""

# ╔═╡ 329d30f1-797e-4522-9c20-e60d35079f5f
adm = read_csv("$(OUTPUTDIR)/tutorial_adm.csv", DataFrame)

# ╔═╡ f550da45-1202-4f9d-9f0b-b96d5c929f58
let
	fig = Figure()
	time = adm[!, 1]

	ax = Axis(fig[1, 1], title="Age-Depth Model", xlabel="t [Myr]", ylabel="depth (m)")

	for i = 1:3
		lines!(ax, time, adm[!, i+1], label="location $i")
	end

	fig[1,2] = Legend(fig, ax)
	fig
end

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
# ╠═ebcb24b3-1ea3-49a8-a6d8-bf1f2cee657e
# ╠═3c4cef70-df77-46ba-b623-fd46b5500e51
# ╟─0ce8de55-3304-431d-a2aa-110b46a25c9b
# ╟─9babc2b0-9c26-11ef-3459-0d113ec3c402
# ╟─17722d8b-baca-4f16-981f-1501c734a95f
# ╟─22ec7e16-b8c0-414d-9700-52bf379e1051
# ╟─3b7fef8b-efb9-467d-b6db-f7cfa132be69
# ╠═bcea7127-3c21-4c35-af42-3d2c71464409
# ╟─cf61bc3f-a20a-45b7-a885-22b70075fc42
# ╟─68fac1d8-f402-429e-90a4-25fcfa188c2e
# ╠═d72d7e42-8392-44a0-a8b3-d59475be8dc7
# ╠═325e3d04-2ff2-4c27-91bf-265820ac6763
# ╟─9aafba01-fc4c-4dc1-85b6-9f33a4cfc77a
# ╠═b3b271cb-143f-44ba-a593-80b9e6c96392
# ╟─1d5cf6cc-745d-4a5a-80ae-b1b6c057af0b
# ╟─c974cb9a-e1c0-4402-880c-7990d217da89
# ╟─316b049d-698d-4d5e-9c18-73701ef8b492
# ╠═74f4674f-dbea-44ad-8d54-9861b35139cd
# ╟─66e30189-ae72-4ec1-b7bd-1136ddfce2ee
# ╠═4fe0f485-2db0-4685-b5b9-e9ba6009e4a6
# ╠═5432b762-50fa-4ac4-97e6-0477236cb94a
# ╠═e118a117-9a00-4589-906d-c31d2057bcef
# ╟─8f883ea5-d90e-41e7-9809-3f170183a640
# ╠═51f2b3db-7b2d-4d62-8c3a-6608c01bffb7
# ╠═b6ae4907-ac04-4a37-975f-42d85d572f91
# ╠═342490a7-685f-47f4-9889-a56927713f85
# ╟─5e914211-fd1f-45e0-9ecd-442d819684a2
# ╟─a39a9476-6a7f-435a-b882-36ecf618369e
# ╟─620627e8-4bf6-400d-ab6f-9a7a9a6ec040
# ╠═37cacf0c-c33a-4a11-b47e-fe2dc24d4eff
# ╠═90fdaca0-2efa-4230-8af3-177a6d94639c
# ╠═281df21c-aa0a-41fb-a20d-45f7b4eb3532
# ╟─a2689006-cb86-4956-999a-7adea546abdd
# ╟─a81dba9b-8f12-4b72-9ea8-993d6c5e501b
# ╠═7d72a3bb-5ed8-4ebf-a45d-59e1bb267e9a
# ╠═30245436-9731-4cae-a150-fcba4360527e
# ╟─f6776dbe-8dc5-4c01-a034-cc47d0c34e8a
# ╟─5c90b464-858a-4a5d-a2c7-fda775057ef6
# ╠═68f2dd79-6777-4267-ae67-5167f764b7b9
# ╟─2e168646-5b9e-437b-b0b2-d637a4beb577
# ╠═d0bfae1a-deef-4625-a242-0a9b899bf83d
# ╠═68ab2b1c-67fb-48f1-ac81-c7efac04d96a
# ╟─43405e49-2739-4e79-8204-e489db6c1fd5
# ╠═ae3a74ea-f6b5-442c-973b-1cec48627968
# ╠═8946d268-4407-4fe4-86ae-67b3a37b34be
# ╠═0fd7ea33-ff06-4355-9278-125c8ed66df4
# ╟─a3e5f420-a59d-4725-8f8f-e5b8f06987db
# ╠═3723832f-344c-4471-acb0-cef7d4e5ca94
# ╠═ec4f0021-12af-4f8c-8acb-970e6820d2a4
# ╠═7f05b817-f933-4280-b2ed-ae318a535123
# ╟─2a24237e-5c1f-47e1-8f33-cca3ef563930
# ╠═31e7c759-f980-4618-be90-892865751e58
# ╠═61ae751b-0c05-4e96-8be9-3f85cb6afc51
# ╠═329d30f1-797e-4522-9c20-e60d35079f5f
# ╠═f550da45-1202-4f9d-9f0b-b96d5c929f58
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
