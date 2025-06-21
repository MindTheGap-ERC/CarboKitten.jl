### A Pluto.jl notebook ###
# v0.20.10

#> [frontmatter]
#> title = "CarboKitten Tutorial"
#> date = "2024-11-18"
#> tags = ["geoscience", "stratigraphy"]
#> description = "a first tutorial on working with CarboKitten and Julia"
#> 
#>     [[frontmatter.author]]
#>     name = "Johan Hidding"
#>     [[frontmatter.author]]
#>     name = "Xianyi Liu"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ abfcbd95-3323-4469-a568-05565675613e
using Pkg; Pkg.activate("../../workenv")

# ╔═╡ bcea7127-3c21-4c35-af42-3d2c71464409
using CarboKitten

# ╔═╡ 325e3d04-2ff2-4c27-91bf-265820ac6763
using CarboKitten.Models: ALCAP

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

# ╔═╡ 94aea87e-89e0-4aa3-816e-397873fffc77
using CarboKitten.Visualization: production_curve

# ╔═╡ 3723832f-344c-4471-acb0-cef7d4e5ca94
using CarboKitten.Export: data_export, CSV

# ╔═╡ 31e7c759-f980-4618-be90-892865751e58
using DataFrames

# ╔═╡ 61ae751b-0c05-4e96-8be9-3f85cb6afc51
using CSV: read as read_csv

# ╔═╡ dec0cd85-bbd7-4a74-83a8-b9425e053f86
using CarboKitten.Export: read_column

# ╔═╡ e892bc37-81d3-4b8f-9486-de0d717cd67f
using CarboKitten.Visualization: stratigraphic_column!

# ╔═╡ 03d95080-35ce-4708-8306-dd5f8dc8c3c0
using Interpolations

# ╔═╡ 9a044d3d-fb41-4ffe-a3ad-5acd94aa6ac6
using DelimitedFiles

# ╔═╡ 950f8d1e-59d8-4485-91fe-2baa21478833
using RCall

# ╔═╡ 17501c93-f432-4f1a-b815-5ac9c5a29f8f
using CarboKitten.DataSets: miller_2020

# ╔═╡ 99d20556-eefb-4597-88a3-0b61bd3cb5c8
using CarboKitten.Boxes: axes as box_axes

# ╔═╡ e1c9135e-4906-4bd1-ba9a-058f7de7d5ac
using CarboKitten.Visualization: sediment_profile!

# ╔═╡ 5eb8bf63-de51-4f4f-ba23-53e4a68e5762
using CarboKitten.Export: age_depth_model, read_slice, DataColumn

# ╔═╡ 44308ced-efd2-42fd-94ab-baa178352ed9
begin
	using ShortCodes
	
	refs = Dict(
		:pomar2001 => DOI("10.1046/j.0950-091x.2001.00152.x"),
		:church2017 => DOI("10.1002/2016WR019675"),
		:burgess2013 => DOI("10.1016/j.cageo.2011.08.026"),
		:bosscher1992 => DOI("10.1111/j.1365-3091.1992.tb02130.x"),
		:paola1992 => DOI("10.1111/j.1365-2117.1992.tb00145.x"),
		:miller2020 => DOI("10.1594/PANGAEA.923126"))

	for d in values(refs)
		ShortCodes.getdoi(d)
	end

	function cite(ref)
		r = refs[ref]
		authors = [strip(split(a, ",")[1]) for a in split(r.author, ";")]
		author = if length(authors) == 1
			authors[1]
		elseif length(authors) == 2
			"$(authors[1]) & $(authors[2])"
		else
			"$(authors[1]) et al."
		end
		year = r.year
		"[$(author) ($(year))](https://doi.org/$(r.doi))"
	end

	sort(refs |> values |> collect, by=r->[r.year, r.author])
end

# ╔═╡ ebcb24b3-1ea3-49a8-a6d8-bf1f2cee657e
begin
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 0ce8de55-3304-431d-a2aa-110b46a25c9b
md"""
This CarboKitten tutorial is aimed at people that are completely new to Julia. Please follow the button on the top-right to download this notebook, and put it in an new directory.
"""

# ╔═╡ 9babc2b0-9c26-11ef-3459-0d113ec3c402
md"""
# Installing

Please install Julia from the following webpage: [https://julialang.org/downloads/](https://julialang.org/downloads/).

We will use [Pluto](https://plutojl.org/) to do our tutorial. This is a notebook interface (similar to Jupyter) that is easy to use and has a strong focus on reproducibility.

Start the Julia REPL (read-eval-print loop), either from the start menu (on Windows), the Launchpad (on Mac), or in a terminal, by typing `julia` and pressing Enter. You should see a colorful welcome message and a prompt for input:
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
You can install Pluto by typing `using Pluto`, pressing Enter and then answering `y` to the prompted question.

```juliarepl
julia> using Pluto
 │ Package Pluto not found, but a package named Pluto is available from a
 │ registry. 
 │ Install package?
 │   (@v1.11) pkg> add Pluto 
 └ (y/n/o) [y]: 
```

Please type `Pluto.run()` and press Enter to start the Pluto Notebook as follows:

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

After a while you should see a pop-up window.


Open the notebook that you downloaded. 
"""

# ╔═╡ d573309f-c99a-43e9-a89f-083ef4ade5d8
md"""
!!! danger "☕ Coffee time! ☕"
	When you open the notebook and click "Run notebook code" at the top-right corner, Pluto will download and compile all the needed packages. Depending on your machine this may take 5-10 minutes. We would love to hear your feedback: please share it via (this form)[https://forms.office.com/e/8RnBXqxq2D].
"""

# ╔═╡ 3b7fef8b-efb9-467d-b6db-f7cfa132be69
md"""
## Install CarboKitten

In a Pluto notebook, Julia packages are installed by `using` them. For example, we enter `using CarboKitten` here to install CarboKitten.jl.
"""

# ╔═╡ cf61bc3f-a20a-45b7-a885-22b70075fc42
md"""
All packages used and their versions are stored inside the notebooks. When you run the notebook on a different computer, the same packages will be used.
"""

# ╔═╡ 68fac1d8-f402-429e-90a4-25fcfa188c2e
md"# Running an existing model"

# ╔═╡ 002cb6d7-ee29-408f-a289-36ab913c7f85
md"""
## Import the model definition
"""

# ╔═╡ 545a6a8d-70d5-470a-a615-4305efa0ecd1
md"""
This is a built-in example model with default values of parameters. It uses a simple sinusoidal curve to represent sea level.
"""

# ╔═╡ dd4cde67-4329-4135-80d7-1c8950404349
Markdown.parse("""
We have now imported the `ALCAP` model, which stands for **A**ctive **L**ayer, **C**ellular **A**utomaton driven **P**roduction. This model uses the CA proposed by $(cite(:burgess2013)), production model by $(cite(:bosscher1992)) and an Active Layer diffusive transport approach (inspired on $(cite(:paola1992)), see $(cite(:church2017)) for a concise overview).
""")

# ╔═╡ 9aafba01-fc4c-4dc1-85b6-9f33a4cfc77a
md"""
## Set the output directory

Please make sure to set the output directory to a convenient place. If you downloaded this notebook to an empty directory, using `"."` would be a good choice.
"""

# ╔═╡ b3b271cb-143f-44ba-a593-80b9e6c96392
OUTPUTDIR = "../../data/output"

# ╔═╡ 316b049d-698d-4d5e-9c18-73701ef8b492
md"""
!!! tip "Disabled cells"
    The cell below is disabled because it will take a minute or so to run. Click the cell-menu at the top right of the cell to enable it.
"""

# ╔═╡ 687e5908-03d7-44ee-a9d8-83f32cf9e447
md"""The `run_model` command runs a model with given input and writes the output to an HDF5 file."""

# ╔═╡ 74f4674f-dbea-44ad-8d54-9861b35139cd
# example_output = run_model(Model{ALCAP}, ALCAP.Example.INPUT, "$(OUTPUTDIR)/example.h5")
example_output = "$(OUTPUTDIR)/alcap_example.h5"

# ╔═╡ 19da029b-8752-4177-8ba4-cc2097adec95
md"""
## Plotting a cross-section
CarboKitten has a set of built-in visualizations. Here we plot a cross-section of the platform along the onshore-offshore gradient.
"""

# ╔═╡ 66e30189-ae72-4ec1-b7bd-1136ddfce2ee
md"""
!!! tip "Makie"
	We use [Makie](https://makie.org) to do our plotting. If you'd like to learn more about Makie, check their extensive documentation.
"""

# ╔═╡ e118a117-9a00-4589-906d-c31d2057bcef
summary_plot(example_output)

# ╔═╡ 56765b03-9d25-49c6-9aec-75e1e32e6a43
md"""
The first sub-diagram on the top left shows the cross-section of the simulated carbonate-platform. Moving clockwise, the second diagram shows a topographic overview of the entire simulation in 3D. The third diagram show the production rate used in this simulation. The fourth diagram shows the sea-level curve. On the bottom left there are two chronostratigraphic (Wheeler's) diagrams, one for sedimentation rate, the other showing the dominant facies type.
"""

# ╔═╡ c974cb9a-e1c0-4402-880c-7990d217da89
md"""
!!! info "Exercise: String interpolation"
	Notice the `"$(...)"` syntax: this is called string interpolation: the contents of the `OUTPUTDIR` variable are formatted into the resulting string. Try it out with a few examples (inside a `let` block so the variables don't leak):

	```julia
	let
		name = "Peter"
		"Hello, $(name)!"
	end
	```

	What happens if you change `name` into a number?
"""

# ╔═╡ 8f883ea5-d90e-41e7-9809-3f170183a640
md"""
# Creating your own model

In our second example we walk through the entire configuration. We start by defining the shape of our box. Overall, we need to define the grids (simulation box), the time (simulation time), the sea-level curve, the initial topography, the facies and the other to start the simulation.

## Simulation box
"""

# ╔═╡ 342490a7-685f-47f4-9889-a56927713f85
box = Box{Coast}(grid_size = (50, 50), phys_scale = 300.0u"m")

# ╔═╡ 5e914211-fd1f-45e0-9ecd-442d819684a2
md"""
There is a lot to unpack here: `Coast` indicates a box topology that is periodic in the y-direction and mirrored in the x-direction, 
which is ideal for simulating a strip of coastline. Other choices could be `Periodic{2}` or `Reflected{2}`. 
The `grid_size` argument is the size of our box in pixels, and `phys_scale` is the size of one pixel with units of length.
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
md"""Notice that the Δt property was automatically converted to Myr, which is the unit that we use internally.
"""

# ╔═╡ a81dba9b-8f12-4b72-9ea8-993d6c5e501b
md"""
## Sea level

In this example, the sea level is given as a function. The function should take quantities in units of time and return a quantity in units of length.

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
	a1 = 4.0u"m"
	p1 = 0.2u"Myr"

	return a1 * sin(2π * t / p1)
end

# ╔═╡ 30245436-9731-4cae-a150-fcba4360527e
let
	t = time_axis(time)   # obtain a vector of times
	sl = sea_level.(t)	    # get sea_level for all those times
	fig, ax = lines(t |> in_units_of(u"Myr"), sl |> in_units_of(u"m"))
	ax.xlabel = "time [Myr]"
	ax.ylabel = "sea level [m]"
	fig
end

# ╔═╡ f6776dbe-8dc5-4c01-a034-cc47d0c34e8a
md"""
!!! info "Exercise: expand the `sea_level` function"
	Try to add a second mode of oscillation into the sea level curve. Choose your own amplitude and period.
"""

# ╔═╡ 5c90b464-858a-4a5d-a2c7-fda775057ef6
md"""
## Initial topography
In CarboKitten, the initial topography is a function of both $x$ and $y$ coordinates, giving the initial sea floor height for each point. Here we specify a simple ramp in the $x$ direction.
"""

# ╔═╡ 68f2dd79-6777-4267-ae67-5167f764b7b9
function initial_topography(x, y)
	return -x / 300.0
end

# ╔═╡ 2e168646-5b9e-437b-b0b2-d637a4beb577
md"""
## Facies definitions

Three biological facies are distinguished based on sediment produced by three carbonate factories:

- T(Tropical) factory (max growth rate = 500.0 m/Myr)
- M(Mounds) factory (max growth rate = 400.0 m/Myr)
- C(Cool water) factory (max growth rate = 100.0 m/Myr)

Their growth rates depend on the water depth, as parametrized by the Bosscher & Schlager 1992 model,

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right),$$

where $g_m$ is the maximum growth rate, $I_0$ is the insolation, $I_k$ the saturation intensity, $k$ the extinction coefficient, and $w$ the water depth.

Because our model also includes diffusive transport, we also specify a diffusion coefficient for each facies type. The units are a bit weird, but the intuition should be: larger values for smaller particle size.

When we inpsect the displayed `Facies` values, we see that also the default viability and activation ranges are set to govern the cellular automaton, following Burgess 2013.
"""

# ╔═╡ 68ab2b1c-67fb-48f1-ac81-c7efac04d96a
const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

# ╔═╡ 43405e49-2739-4e79-8204-e489db6c1fd5
md"""
## Other parameters

The subsidence rate is set to a constant 50 m/Myr and insolation to 400 $W/m^2$. The other parameters are related to the transport model.
"""

# ╔═╡ ae3a74ea-f6b5-442c-973b-1cec48627968
input = ALCAP.Input(
    tag="tutorial",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> sea_level(t),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)


# ╔═╡ 549c852c-44d8-483d-b54c-cdf68ff6020b
production_curve(input)

# ╔═╡ 8946d268-4407-4fe4-86ae-67b3a37b34be
tutorial_output = run_model(Model{ALCAP}, input, "$(OUTPUTDIR)/tutorial.h5")

# ╔═╡ 0fd7ea33-ff06-4355-9278-125c8ed66df4
summary_plot(tutorial_output)

# ╔═╡ 39e681ae-3ce4-41eb-a48a-9a42a4667183
md"""
!!! info "Exercise: sea level again"
	With the `run_model` instruction above still active, modify your input sea level curve again, like you did in the previous exercise. Re-run the `summary_plot` command again to see the result.
"""

# ╔═╡ a3e5f420-a59d-4725-8f8f-e5b8f06987db
md"""
# Extracting model output into CSV files

HDF5 is a very versatile data format, and is readable from every major computational platform (e.g. MatLab, Python, R). However, sometimes you may want to process your output data further using existing tools that read CSV data.
"""

# ╔═╡ ec4f0021-12af-4f8c-8acb-970e6820d2a4
export_locations = [(11, 25), (26, 25), (41, 25)]

# ╔═╡ e80aaf02-c5bf-4555-a82d-b09dcf785381
md"""
In this case, it exports the 11th, 26th and 41th grid in the direction towards to the deep sea. That is to say, the 11th is proximal to the land while the 41th is distal. Given that each grid is with fixed size of 300 m × 300 m, this indicates that we are extracting information from locations at 3 km, 7.5 km and 12 km away from the land. 
"""

# ╔═╡ 7f05b817-f933-4280-b2ed-ae318a535123
# ╠═╡ disabled = true
#=╠═╡
data_export(CSV(export_locations,
	:sediment_accumulation_curve => "$(OUTPUTDIR)/tutorial_sac.csv",
	:age_depth_model => "$(OUTPUTDIR)/tutorial_adm.csv"),
	tutorial_output)
  ╠═╡ =#

# ╔═╡ 2a24237e-5c1f-47e1-8f33-cca3ef563930
md"""
!!! info "Exercise: Meta data"
	Next to `:sediment_accumulation_curve` and `:age_depth_model`, we have implemented functions to extract `:stratigraphic_column` and `:metadata`.

	The `:metadata` target is special, since it doesn't write to CSV but to TOML. Add the `:stratigraphic_column` and `:metadata` targets to the export a list and inspect the result.
"""

# ╔═╡ adf67033-5941-4be6-bb17-c0958348b905
md"""
## Plot age-depth models
"""

# ╔═╡ e26defa5-10ff-4275-bae9-768f7fb8d9ba
md"""
We use the `DataFrames` package to process the CSV file, and again use Makie to visualise the data.
"""

# ╔═╡ e17f35e7-8d09-4da1-880f-563bc49b364c
md"""
Using the following command to read the csv file:
"""

# ╔═╡ 329d30f1-797e-4522-9c20-e60d35079f5f
# ╠═╡ disabled = true
#=╠═╡
adm = read_csv("$(OUTPUTDIR)/tutorial_adm.csv", DataFrame)
  ╠═╡ =#

# ╔═╡ add6a25b-d948-4cd0-9412-56752793ca4b
md"""
And plot the Age-Depth Model:
"""

# ╔═╡ f550da45-1202-4f9d-9f0b-b96d5c929f58
#=╠═╡
let
	fig = Figure()
	time = adm[!, 1]

	(xaxis, _) = box_axes(input.box)
	ax = Axis(fig[1, 1], title="Age-Depth Model", xlabel="t [Myr]", ylabel="depth (m)")
	

	for i = 1:3
		xpos = uconvert(u"km", xaxis[export_locations[i][1]])
		lines!(ax, time, adm[!, i+1], label="$(xpos) offshore")
	end

	fig[1,2] = Legend(fig, ax)
	fig
end
  ╠═╡ =#

# ╔═╡ 64c3ce44-de95-4f7b-954d-52f743fc5033
md"""
## Plot the stratigraphic column
"""

# ╔═╡ 0a5da056-ca6a-4d90-b2a4-fb84ae5b7da2
let
	fig = Figure()

	y = export_locations[1][2]
	(xaxis, _) = box_axes(input.box)
	
	for i = 1:3
		header, column = read_column(tutorial_output, export_locations[i]...)

		xpos = uconvert(u"km", xaxis[export_locations[i][1]])
		ax = Axis(fig[1,i], title="$(xpos) offshore", ylabel="height (m)", xticksvisible=false, xtickformat="", width=70)
		stratigraphic_column!(ax, header, column)
	end
	
	fig
end

# ╔═╡ d7a4eb34-6b9e-4779-900d-493a853a6199
md"""
# Try a sea-level curve that is modulated by orbital cycles

You may need to import an external source, for example, a txt, a excel, or a csv file. Herein, we use a txt file, and therefore we need to add dependencies of 'DelimitedFiles'.

If you did not prepare your sea-level curve, we provide one and you could download it from emails. You could download it to your own computer, and replace the SL_PTH with the real pathway that you stored the txt file. It is recommended store the notebook and your sealevel curve at the same directory.

"""

# ╔═╡ 7f6c8d9e-db15-48fc-a631-662d56f87197
SL_PTH = "sl.txt"

# ╔═╡ 41e84b67-f5f6-4b24-943d-b0bfe50ca2bb
begin
		input_sl = DataFrame(readdlm(SL_PTH, '\t', header=false),:auto)
		input_sl_clean = DataFrame(filter(row -> all(x -> string(x) != "", row), eachrow(input_sl)))
		sea_level_ob = input_sl_clean[:,2] .* 1.0u"m"
		time_ob = input_sl_clean[:,1].* 1.0u"Myr"
		Interpolated_SL = linear_interpolation(time_ob,sea_level_ob)
end

# ╔═╡ 5043ab85-a7ff-4283-9d71-071c470f4d6e
input_ob = ALCAP.Input(
    tag="tutorial",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=3900,
        write_interval=1),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level= Interpolated_SL,
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

# ╔═╡ bd208d68-f565-4d05-89f7-f2e572a87f04
# ╠═╡ disabled = true
#=╠═╡
own_sealevel_output = run_model(Model{ALCAP}, input_ob, "$(OUTPUTDIR)/tutorial_ownsealevel.h5")
  ╠═╡ =#

# ╔═╡ b2286eff-a49e-4596-b68c-1a5363c9477f
# ╠═╡ disabled = true
#=╠═╡
data_export(CSV(export_locations,
	:sediment_accumulation_curve => "$(OUTPUTDIR)/tutorial_ownsealevel_sac.csv",
	:age_depth_model => "$(OUTPUTDIR)/tutorial_ownsealevel_adm.csv"),
	own_sealevel_output)
  ╠═╡ =#

# ╔═╡ 7f7dd756-d2e4-4eeb-8364-ea750a0aecc3
md"""
We assume we can find sufficient forams (i.e., ignore sampling bias) from the pseudo strata. We then 'crush' these 'forams' to obtain their Oxygen isotope values, trying to understand how the temperature changes through time. In this case, what kind of oxygen isotope curve did you expect?

First we constructed a pseudo Oxygen isotope curve based on the sea-level curve (e.g.,assuming low-sea-level is glaciation and high O isotope composition). We do not consider the offset of phases.

We then use the ADM extracted from the carbonate platform to test how much of the original signals would be preserved. Could the preserved signals still refect all of the orbital cycles? 

"""

# ╔═╡ 4b83e9e9-6404-4e52-8662-7ec56ccc7fa6
let
	O_magnitude = 4 #arbitrary
	sl_magnitude = (max(sea_level_ob...) - min(sea_level_ob...)) ./ u"m"
	d18O_original = sea_level_ob ./ sl_magnitude .* O_magnitude
	
	function determine_preservation(ADM::Vector)
		ADM_diff = diff(ADM)
		preservation_potential = Float64[]
			for i in 1:length(ADM_diff)
				if ADM_diff[i] > 0.0001 * u"m" #arbitrary
					push!(preservation_potential, 1.0)
				else
					push!(preservation_potential, NaN)
				end
			end
		return preservation_potential
	end

	function d18O_preserve(preservation_potential::Vector,d18O::Vector,time_input::Vector, time_output::Vector)
		d18O_inter = linear_interpolation(time_input,d18O)
		d18O_out = d18O_inter.(time_output)
		return preservation_potential .* d18O_out[1:end-1]
	end
	
	adm = read_csv("$(OUTPUTDIR)/tutorial_ownsealevel_adm.csv", DataFrame)
	adm_shallow = adm[:,2] .* u"m"
	adm_deep = adm[:,end] .* u"m"
	time = adm[:,1].* u"Myr"

	d180_shallow = d18O_preserve(determine_preservation(adm_shallow), d18O_original,time_ob,time)
	d180_deep = d18O_preserve(determine_preservation(adm_deep), d18O_original,time_ob,time)

	writedlm("$(OUTPUTDIR)/d180_shallow.txt", d180_shallow, '\t')
	writedlm("$(OUTPUTDIR)/d180_deep.txt", d180_deep, '\t')

end

# ╔═╡ bcd404d1-d2a8-4b8c-8fe8-16abb9215442
md"""
Similarly, shall we also try to change the insolation from a const (400 W/m2) to a real insolation curve, and see how it works?
"""

# ╔═╡ b798b991-5af9-4c4c-b795-88e8087c6c4d
md"""
We first import a curve from R package palisol
"""

# ╔═╡ e56a4012-d80a-4535-a575-ad80c9d28610


# ╔═╡ eb3908ca-2d33-437d-9a19-a65efa20da3e
R"""
if (!require("palisol")) {
    install.packages("palisol", repos = "https://cran.r-project.org")
}
"""


# ╔═╡ b1fee3c8-bd01-4f07-b2ed-18730999967f
R"""
library(palisol)
time_start <- 8e5  
time_end <- 0      
time_step <- 2e2  
times <- seq(time_end, time_start, time_step)

param_la04 = astro(times, solution = la04,  degree = TRUE)

orbit <- list()
insolation <- list()
lat_degree = 25

for (t in 1:length(times)) {
  orbit[[t]] <- list(
    eps = param_la04[1] * pi / 180, 
    ecc = param_la04[2], 
    varpi = (param_la04[t + 2] - 180) * pi / 180
  )
  
  insolation[[t]] <- Insol(
    orbit[[t]], 
    long = pi / 2, 
    lat = lat_degree * pi / 180, 
    S0 = 1365, 
    H = NULL
  )
}

insolation = inso_values <- unlist(insolation)

"""


# ╔═╡ 9ad42bc3-0f74-48fa-8668-1d7a9555d992
begin
	annual_insol_julia = rcopy(R"insolation")
	function get_inso()
		interpolated_inso = linear_interpolation(times,insolation)
	end
end

# ╔═╡ af878041-b78b-4e62-8aed-826506935f89
# ╠═╡ disabled = true
#=╠═╡
input_inso = ALCAP.Input(
    tag="tutorial",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=3900,
        write_interval=1),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level= Interpolated_SL,
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation= get_inso() * 1.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)
  ╠═╡ =#

# ╔═╡ 459feada-ad61-4939-b733-30aa007bb026
# ╠═╡ disabled = true
#=╠═╡
own_inso_output = run_model(Model{ALCAP}, input_inso, "$(OUTPUTDIR)/tutorial_ownsinso.h5")
  ╠═╡ =#

# ╔═╡ 7c100b7f-7661-4250-b999-fdd3f32bf80b
# ╠═╡ disabled = true
#=╠═╡
data_export(CSV(export_locations,
	:sediment_accumulation_curve => "$(OUTPUTDIR)/tutorial_owninso_sac.csv",
	:age_depth_model => "$(OUTPUTDIR)/tutorial_owninso_adm.csv"),
	own_inso_output)
  ╠═╡ =#

# ╔═╡ c2166805-da62-4adf-8514-fd28924d115e
Markdown.parse("""
# Using an empirical sea-level curve

We've already seen how we can read CSV files. To make this part a bit easier, CarboKitten ships with the Cenozoic sea level dataset by $(cite(:miller2020)).
""")

# ╔═╡ 71f2c3ad-80ea-4678-87cf-bb95ef5b57ff
miller_df = miller_2020()

# ╔═╡ 04824ca5-3e40-40e7-8211-990c18af21e6
let
	fig = Figure(size=(1000,300))
	ax = Axis(fig[1, 1]; xlabel="time (Ma BP)", ylabel="sealevel (m)")
	
	for ref in levels(miller_df.reference)
		subset = miller_df[miller_df.reference .== ref,:]
		lines!(ax, subset.time |> in_units_of(u"Myr"), subset.sealevel |> in_units_of(u"m"), label=ref)
	end
	fig[1, 2] = Legend(fig, ax)
	fig
end

# ╔═╡ 31a137be-f726-47f9-be2d-a6f91e4084dc
md"""
We can use the `refkey` field to select on a specific dataset. Inspect the unique values using the `levels` function.
"""

# ╔═╡ 3d550592-93a0-4eeb-891e-f690a3b1d686
levels(miller_df.refkey)

# ╔═╡ c0970bfe-db0a-46cc-a3e2-13e1ee78ee8e
md"""
!!! info "Exercise: data selection"
	Read the above code to plot the data in Miller et al. 2020, focussing on creating a subset of the data. Use the `refkey` column to select only the data from Lisiecki et al. 2005, and store it in `lisiecki_df`.

	Plot the result.
"""

# ╔═╡ 6f60233c-739d-40d5-8c43-d30647006351
md"""
As the observed data is not regularly sampled, we need to interpolate the dataset. The `Interpolations` package has a convenient linear interpolation function. It does require input values to be sorted.
"""

# ╔═╡ e2562586-d03a-4b6e-9ef6-aad012f2be9f
# ╠═╡ disabled = true
#=╠═╡
using Interpolations
  ╠═╡ =#

# ╔═╡ 54d26634-765f-4291-8a82-34e2e6e2fa09
sort!(miller_df, [:time]);

# ╔═╡ 1977eab2-f277-41c4-ac00-826bafa8de4d
miller_sea_level = linear_interpolation(miller_df.time, miller_df.sealevel);

# ╔═╡ ca0498b0-c425-4949-b1d9-03df4067db2e
md"""
!!! info "Exercise: run the model"
	Create a sea level function by interpolating the Lisiecki et al. data set. Make sure to start the simulation at t=-2.0 Myr, for example:

	```julia
	time = TimeProperties(
		t0 = -2.0u"Myr",
		Δt = 500.0u"yr",
		steps = 2000
	)
	```

	Try to plot the cross-section and adm from this simulation.
"""

# ╔═╡ 3afb006d-b6ff-4ae7-8e60-115d98e9d562
md"""
!!! info "Exercise: Input your own sea-level curve"
	We've seen how to load CSV data, and how we can interpolate tabular data to create a new sea-level curve. Can you run the model using your own sea-level data?
"""

# ╔═╡ 04a28420-a758-4597-a679-b675844d9c7e
md"""
# Interactive Exploration
"""

# ╔═╡ cfb15723-694b-4dc1-bd37-21d17074ab98
md"""
In Pluto we can create interactive visualizations. Sliding `y` will cause a reload of the data slices, which can be a bit slow.

!!! tip "Macros"
	In Julia, the `@` character indicates the use of a macro. It's Ok not to understand: macro is a different word for magic smoke screen.

	If you'd like to learn more about interactivity in notebooks, check out the [`PlutoUI` demo notebook](https://featured.plutojl.org/basic/plutoui.jl).
"""

# ╔═╡ 1fc8dffe-8e56-4da9-a7e0-07793d1d0455
@bind y PlutoUI.Slider(1:50, default=25)

# ╔═╡ 8dfe7d41-2915-48b7-9866-b36ae7fe6443
y

# ╔═╡ 8e0bdea3-77a7-452f-a3ca-b2a0b5bbf5f0
md"""
We load the data using the `read_slice` function from `CarboKitten.Export`. This extracts a slice along the $x$-axis (onshore-offshore profile).
"""

# ╔═╡ a3485ec5-bdf1-4577-9d31-3ea21eba6a53
header, data_slice = read_slice(tutorial_output, :, y);

# ╔═╡ fbf6306d-de6c-4581-9f19-4ac066162709
md"""
Within this slice we can select a single stratigraphic column by setting the $x$ coordinate."""

# ╔═╡ 1d952bc7-4375-4f6b-a47a-a2b0ea915360
@bind x PlutoUI.Slider(1:50)

# ╔═╡ eb4a1b28-4ba4-4cd1-9490-e27105b19921
x

# ╔═╡ c820f883-fa5b-4a43-a03f-e6961dbe4001
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	(_, yaxis) = box_axes(input.box)
	ypos = yaxis[y]
	sediment_profile!(ax, header, data_slice)
	ax.title = "sediment profile @ y = $(uconvert(u"km", ypos))"
	vlines!(ax, header.axes.x[x] |> in_units_of(u"km"), color=:black)
	fig
end

# ╔═╡ 5ca6c8ed-cfd1-478d-a9fe-67fb4ff2ef0f
let
	fig = Figure()
	(xaxis, yaxis) = box_axes(input.box)
	xpos = uconvert(u"km", xaxis[x])
	ypos = uconvert(u"km", yaxis[y])
	ax = Axis(fig[1, 1], title="age depth model @ (x = $(xpos) offshore, y = $(ypos))", ylabel="height [m]", xlabel="time [Myr]")
	column = DataColumn(
		(x, y),
		data_slice.disintegration[:, x, :],
        data_slice.production[:, x, :],
        data_slice.deposition[:, x, :],
        data_slice.sediment_elevation[x, :])
	adm = age_depth_model(column.sediment_elevation) |> in_units_of(u"m")
	lines!(ax, header.axes.t |> in_units_of(u"Myr"), adm)

	ax2 = Axis(fig[1, 2], title="strat. col.", width=70, xtickformat="", xticksvisible=false)
	stratigraphic_column!(ax2, header, column)
	fig
end

# ╔═╡ c85294f2-3508-48c0-b67d-d82df646ae33
Markdown.parse("""
# Task: Platform Morphology

Researchers have found different morphologies of carbonate platforms. For example, some of them show a 'steep cliff' (rimmed-shelf) while some others show a 'smooth' ramp. See for instance $(cite(:pomar2001)).

Different carbonate producers (i.e., T, M, C) produce carbonate with different production rate under different water-depth. Could the production rate be a key controller for the morphology of the carbonate platform? That is to say, can you vary the parameters as illustrated in section on Facies Definitions to try to produce two carbonate platforms, one with 'rimmed-shelf' and another one with ramp'?
""")

# ╔═╡ 754f7cf0-f3a7-4cff-b20e-ed16874860c7
md"""
# Bibliography
"""

# ╔═╡ Cell order:
# ╠═abfcbd95-3323-4469-a568-05565675613e
# ╟─0ce8de55-3304-431d-a2aa-110b46a25c9b
# ╟─9babc2b0-9c26-11ef-3459-0d113ec3c402
# ╟─17722d8b-baca-4f16-981f-1501c734a95f
# ╟─22ec7e16-b8c0-414d-9700-52bf379e1051
# ╟─d573309f-c99a-43e9-a89f-083ef4ade5d8
# ╟─3b7fef8b-efb9-467d-b6db-f7cfa132be69
# ╠═bcea7127-3c21-4c35-af42-3d2c71464409
# ╟─cf61bc3f-a20a-45b7-a885-22b70075fc42
# ╟─68fac1d8-f402-429e-90a4-25fcfa188c2e
# ╟─002cb6d7-ee29-408f-a289-36ab913c7f85
# ╟─545a6a8d-70d5-470a-a615-4305efa0ecd1
# ╠═325e3d04-2ff2-4c27-91bf-265820ac6763
# ╟─dd4cde67-4329-4135-80d7-1c8950404349
# ╟─9aafba01-fc4c-4dc1-85b6-9f33a4cfc77a
# ╠═b3b271cb-143f-44ba-a593-80b9e6c96392
# ╟─316b049d-698d-4d5e-9c18-73701ef8b492
# ╟─687e5908-03d7-44ee-a9d8-83f32cf9e447
# ╠═74f4674f-dbea-44ad-8d54-9861b35139cd
# ╟─19da029b-8752-4177-8ba4-cc2097adec95
# ╟─66e30189-ae72-4ec1-b7bd-1136ddfce2ee
# ╠═4fe0f485-2db0-4685-b5b9-e9ba6009e4a6
# ╠═5432b762-50fa-4ac4-97e6-0477236cb94a
# ╠═e118a117-9a00-4589-906d-c31d2057bcef
# ╟─56765b03-9d25-49c6-9aec-75e1e32e6a43
# ╟─c974cb9a-e1c0-4402-880c-7990d217da89
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
# ╠═68ab2b1c-67fb-48f1-ac81-c7efac04d96a
# ╠═94aea87e-89e0-4aa3-816e-397873fffc77
# ╠═549c852c-44d8-483d-b54c-cdf68ff6020b
# ╟─43405e49-2739-4e79-8204-e489db6c1fd5
# ╠═ae3a74ea-f6b5-442c-973b-1cec48627968
# ╠═8946d268-4407-4fe4-86ae-67b3a37b34be
# ╠═0fd7ea33-ff06-4355-9278-125c8ed66df4
# ╟─39e681ae-3ce4-41eb-a48a-9a42a4667183
# ╟─a3e5f420-a59d-4725-8f8f-e5b8f06987db
# ╠═3723832f-344c-4471-acb0-cef7d4e5ca94
# ╠═ec4f0021-12af-4f8c-8acb-970e6820d2a4
# ╟─e80aaf02-c5bf-4555-a82d-b09dcf785381
# ╠═7f05b817-f933-4280-b2ed-ae318a535123
# ╟─2a24237e-5c1f-47e1-8f33-cca3ef563930
# ╟─adf67033-5941-4be6-bb17-c0958348b905
# ╟─e26defa5-10ff-4275-bae9-768f7fb8d9ba
# ╠═31e7c759-f980-4618-be90-892865751e58
# ╠═61ae751b-0c05-4e96-8be9-3f85cb6afc51
# ╟─e17f35e7-8d09-4da1-880f-563bc49b364c
# ╠═329d30f1-797e-4522-9c20-e60d35079f5f
# ╟─add6a25b-d948-4cd0-9412-56752793ca4b
# ╠═f550da45-1202-4f9d-9f0b-b96d5c929f58
# ╟─64c3ce44-de95-4f7b-954d-52f743fc5033
# ╠═dec0cd85-bbd7-4a74-83a8-b9425e053f86
# ╠═e892bc37-81d3-4b8f-9486-de0d717cd67f
# ╠═0a5da056-ca6a-4d90-b2a4-fb84ae5b7da2
# ╠═d7a4eb34-6b9e-4779-900d-493a853a6199
# ╠═03d95080-35ce-4708-8306-dd5f8dc8c3c0
# ╠═9a044d3d-fb41-4ffe-a3ad-5acd94aa6ac6
# ╠═7f6c8d9e-db15-48fc-a631-662d56f87197
# ╠═41e84b67-f5f6-4b24-943d-b0bfe50ca2bb
# ╠═5043ab85-a7ff-4283-9d71-071c470f4d6e
# ╠═bd208d68-f565-4d05-89f7-f2e572a87f04
# ╠═b2286eff-a49e-4596-b68c-1a5363c9477f
# ╟─7f7dd756-d2e4-4eeb-8364-ea750a0aecc3
# ╠═4b83e9e9-6404-4e52-8662-7ec56ccc7fa6
# ╟─bcd404d1-d2a8-4b8c-8fe8-16abb9215442
# ╟─b798b991-5af9-4c4c-b795-88e8087c6c4d
# ╠═950f8d1e-59d8-4485-91fe-2baa21478833
# ╠═e56a4012-d80a-4535-a575-ad80c9d28610
# ╠═eb3908ca-2d33-437d-9a19-a65efa20da3e
# ╠═b1fee3c8-bd01-4f07-b2ed-18730999967f
# ╠═9ad42bc3-0f74-48fa-8668-1d7a9555d992
# ╠═af878041-b78b-4e62-8aed-826506935f89
# ╠═459feada-ad61-4939-b733-30aa007bb026
# ╠═7c100b7f-7661-4250-b999-fdd3f32bf80b
# ╠═c2166805-da62-4adf-8514-fd28924d115e
# ╠═17501c93-f432-4f1a-b815-5ac9c5a29f8f
# ╠═71f2c3ad-80ea-4678-87cf-bb95ef5b57ff
# ╟─04824ca5-3e40-40e7-8211-990c18af21e6
# ╟─31a137be-f726-47f9-be2d-a6f91e4084dc
# ╠═3d550592-93a0-4eeb-891e-f690a3b1d686
# ╟─c0970bfe-db0a-46cc-a3e2-13e1ee78ee8e
# ╟─6f60233c-739d-40d5-8c43-d30647006351
# ╠═e2562586-d03a-4b6e-9ef6-aad012f2be9f
# ╠═54d26634-765f-4291-8a82-34e2e6e2fa09
# ╠═1977eab2-f277-41c4-ac00-826bafa8de4d
# ╟─ca0498b0-c425-4949-b1d9-03df4067db2e
# ╟─3afb006d-b6ff-4ae7-8e60-115d98e9d562
# ╟─04a28420-a758-4597-a679-b675844d9c7e
# ╠═99d20556-eefb-4597-88a3-0b61bd3cb5c8
# ╠═e1c9135e-4906-4bd1-ba9a-058f7de7d5ac
# ╠═5eb8bf63-de51-4f4f-ba23-53e4a68e5762
# ╟─cfb15723-694b-4dc1-bd37-21d17074ab98
# ╠═1fc8dffe-8e56-4da9-a7e0-07793d1d0455
# ╠═8dfe7d41-2915-48b7-9866-b36ae7fe6443
# ╟─8e0bdea3-77a7-452f-a3ca-b2a0b5bbf5f0
# ╠═a3485ec5-bdf1-4577-9d31-3ea21eba6a53
# ╟─fbf6306d-de6c-4581-9f19-4ac066162709
# ╠═1d952bc7-4375-4f6b-a47a-a2b0ea915360
# ╠═eb4a1b28-4ba4-4cd1-9490-e27105b19921
# ╟─c820f883-fa5b-4a43-a03f-e6961dbe4001
# ╟─5ca6c8ed-cfd1-478d-a9fe-67fb4ff2ef0f
# ╟─c85294f2-3508-48c0-b67d-d82df646ae33
# ╟─754f7cf0-f3a7-4cff-b20e-ed16874860c7
# ╟─44308ced-efd2-42fd-94ab-baa178352ed9
# ╟─ebcb24b3-1ea3-49a8-a6d8-bf1f2cee657e
