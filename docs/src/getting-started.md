# Getting Started with CarboKitten.jl

This guide takes you from zero to running your first simulation and creating visualizations in CarboKitten.jl, covering installation, writing your own simulation, plotting and exporting data for analysis using other tools.

## Prerequisites and installation

### System requirements

CarboKitten.jl requires **Julia ≥ 1.10** and works across Windows, macOS, and Linux. It is useful to have at least 4GB RAM for model simulations and sufficient disk space for HDF5 output files. For visualization, ensure your system supports OpenGL 3.3+.

If Julia isn't installed, visit https://julialang.org/downloads/ and follow the installation instructions. Verify your installation by running `julia` in your terminal.

### Project setup

**Create a dedicated project directory** with version control to ensure reproducible simulations:

```bash
mkdir MyCarboKittenProject
cd MyCarboKittenProject
git init
julia
```

**Set up the Julia environment** by pressing `]` to enter package mode:

```julia
(@v1.11) pkg> activate .
(MyCarboKittenProject) pkg> add CarboKitten
(MyCarboKittenProject) pkg> add GLMakie  
```

This creates a local `Project.toml` file and installs CarboKitten and its dependencies. In most cases, you will create, modify and run CarboKitten in a dedicated *project*. A Julia project is akin to a project in R Studio or a virtual environment in Python. It keeps the dependencies under control and allows others to reproduce your results. Downloading and compilation may take several minutes on first installation. After this initial long wait, the compiled packages will be cashed on your drive and load faster on subsequent runs.

### Verify installation

Test your installation with this verification sequence:

```julia
using CarboKitten
CarboKitten.init()  
run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
```

The simulation should complete without errors and create an `example.h5` file. You can verify it by plotting a slice of the output. The first time you plot, you need to load GLMakie:

```julia
using GLMakie
using CarboKitten.Visualization
summary_plot("example.h5")
```

A plot window should open showing the carbonate platform cross section. The first rendering may be slow due to compilation, but subsequent plots will be much faster.

## Elements of the model setup you can modify

In this document we refer to running the ALCAP model, CarboKitten's flagship model. Other models are stripped of come of ALCAP's key functionalities and are reproduced for historical reasons or testing purposes.

### Core model components

**Box** represents the spatial domain - the grid where water depth and cellular automata rules simulate facies distributions across the onshore-offshore axis *x* and along-shore axis *y*. The onshore direction is arbitrary and reflects only the convention of showing the shore (higher elevation) on the left-hand side. One can remove it completely by setting `initial_topography` to a flat surface or reverse the convention and set the `initial_topography` to rise with *x*. The `phys_scale` parameters defines the absolute size (in units of length) of one grid cell. `Box` is followed by the indication of the [stencil operations](stencils.md) used on the edges of the grid, e.g. `Box{Coast}`.

**TimeProperties** controls the time when the simulation starts (in units of time), the duration of a time step (in units of time) and the number of time steps (unitless). CarboKitten.jl has been mostly tested for time steps of the order of 100-500 y and time steps wildly different from these may require adjusting other parts of the model, particularly related to transport.

## Basic simulation workflow

### Setting up your first simulation

We strongly recommend **creating a script**, i.e. a saved file, for each of your simulations, even if you're just playing around. You can use any text editor and save the script with the `.jl` extension. Most IDEs are suitable for writing and executing Julia scripts, e.g. VS Code/Codium or R Studio. 

CarboKitten can save model output into the memory or into an HDF5 file, from where it can be read again. The **HDF5 format** (`.h5` files) provides hierarchical, self-documenting data storage with metadata. This format is cross-platform compatible and readable by multiple programming languages, including Python, MATLAB, and R. In the following examples we will use some HDF5 files. They are binary and therefore not suitable to be tracked by git, so you can put them into a separate folder `data/output` and into `.gitignore`. The script won't work if the folder does not exist.

```shell
mkdir data
mkdir data/output
echo "data/output/*.h5" >> .gitignore
```
The easiest way to start modelling is to modify the existing script. The following is the default script that you have previously run when testing the setup, used throughout CarboKitten documentation.

``` julia
#| creates: examples/model/alcap/run.jl
```

You can save this script under a new name in your project and change the parameters there.

## Working with simulation outputs

#### Reading and accessing HDF5 data

You can open the model contents using a general package for handling HDF5 datasets, but we recommend using the functions written specifically for handling CarboKitten outputs. 

This will load the entire model into the memory from an existing file:

``` julia
using CarboKitten.Export: read_volume

header, data = read_volume("data/output/alcap-example.h5", :full)
```
Then the model becomes accessible for examination as described in [Output](https://mindthegap-erc.github.io/CarboKitten.jl/dev/memory-writer/).

In a typical situation of setting up a simulation and adjusting the parameters, we might want to take a look at a cross section of an output. In such case it will be much faster to save and access only one slice of the grid. In the `INPUT` struct we would define then:

```julia
output=Dict(
        :profile => OutputSpec(slice=(:, 25), write_interval=1)),
```

which would save only a slice at the coordinate of y = 25, in the previous example corresponding to the middle of the grid. We can access it using:

```julia
using CarboKitten.Export: read_slice

header, data = read_slice("data/output/alcap-example.h5", :profile)
```

which can then be plotted and saved with:

``` julia
using GLMakie
using CarboKitten.Visualization: sediment_profile

fig = sediment_profile(header, data)
save("data/output/alcap_cross_section_y25.png", fig)
```

If we already have a file containing the entire output, we can extract parts of it:

``` julia
using CarboKitten.Export: read_volume

header, data = read_volume("data/output/alcap-example.h5", :full)
```

#### Exporting and reading CSV files

The example script used above defines what parts of the simulation are exported into CSV. 

``` {.julia #alcap-data_export}

```

CSV does not allow storing arrays with more than 2 dimensions nor metadata, so we don't recommend it as the main format. But you may want to export some data e.g. for easier sharing. 

## Creating visualizations and plots

### Basic plotting workflow


### Custom visualization options


### Plot export formats

Multiple output formats for different use cases:

```julia
fig = summary_plot("simulation.h5")

save("platform_visualization.png", fig, px_per_unit=2)  

save("publication_plot.pdf", fig, pt_per_unit=1.0)
save("scalable_figure.svg", fig)  
```

## Project organization 

You can take inspiration from this proven structure:

```
MyCarboKittenProject/
├── data/                    # Input and output datasets
⎪   └── output/              # Output HDF5 files
├── src/                     # Scripts
├── figs/                    # Generated plots and visualizations
├── Project.toml             # Julia dependencies
├── Manifest.toml            # Exact package versions
├── .gitignore               # Files not to be tracked
├── LICENSE                  # License for your code
└── README.md                # Running instructions
```

## Common issues and solutions

**Installation problems** usually involve Julia version compatibility (ensure ≥ 1.10) or GLMakie graphics issues (update drivers, install OpenGL libraries on Linux). 

**Slow first runs** are normal due to compilation - subsequent runs in the same REPL session will be much faster.

**Memory issues** with large simulations may require careful data management or consideration of model resolution. Consider saving outputs as HDF5 rather than writing them into the memory, using the `write_interval` option to rarify what of the model is saved,

**Genuine bugs** are likely to be found and please don't hesitate to file a bug or suspected bug report in the [issue tracker](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues).

**Incomplete documentation** in spite of best efforts, CarboKitten documentation is not finished. If you're missing something, it is almost certainly our fault. Please don't hesitate to file this in the [issue tracker](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues) or contact one of the developers by email.