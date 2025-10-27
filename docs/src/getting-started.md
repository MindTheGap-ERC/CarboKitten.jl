# Getting Started

## Prerequisites

CarboKitten requires Julia &ge; 1.10. If you haven't installed Julia yet, please follow the [download and install instructions at the Julia homepage](https://julialang.org/downloads/). The recommended way to install Julia is using `juliaup` at [github.com/JuliaLang/juliaup](https://github.com/JuliaLang/juliaup).

## Setting Up a New Project

The instructions below assume you are working in a command line interface (CLI). On Windows, this can be done using PowerShell. Most steps can be performed using a graphical user interface instead, up to the point when you've opened Julia REPL.

### Create a Project Directory

The example here uses `MyCarboKittenProject` as a dummy name. Please replace it with something more meaningful when running these commands on your computer.

```bash
mkdir MyCarboKittenProject
cd MyCarboKittenProject
```

### Initialize Julia Environment

Start Julia:

```bash
julia
```
Create a new project environment
In the Julia REPL, press `]` to enter package mode, then:

```julia
(@v1.11) pkg> activate .
```

You may see a different number, depending on the Julia version you've installed.

After executing it, the prompt should change to the name of your project:

```julia
(MyCarboKittenProject) pkg>
```
In your project folder, you should now see two new files: `Manifest.toml` and `Project.toml`. They contain information about the dependencies you are using and allow reproducing your environment. These files should be tracked by the version control system.

## Installing CarboKitten

Please note that you always install CarboKitten to a given project. This is different from a global installation on your computer. If you want to use CarboKitten in a given project, you always have to do this step. If you have CarboKitten installed in a different project, but not in the current one, you won't be able to use it until you do the installation step:

### Option 1: Install from JuliaHub (Recommended)

From package mode in the Julia REPL:

```julia
(MyCarboKittenProject) pkg> add CarboKitten
```

This will download and install CarboKitten's latest release and its dependencies. 

### Option 2: Install from GitHub

To install the latest development version (i.e. updates since the last release) directly from GitHub:

```julia
(MyCarboKittenProject) pkg> add https://github.com/MindTheGap-ERC/CarboKitten.jl
```

## Installing Unitful

The package allows you to handle units of physical entities and convert between them. It can be installed the same way as option 1 for CarboKitten, from package mode in the Julia REPL:

```julia
(MyCarboKittenProject) pkg> add Unitful
```

Press backspace to exit package mode and return to the Julia REPL.

## Running Your First Model

The steps below propose having two separate scripts for running a model and for plotting its results. This is because plotting requires graphical dependencies, which take additional time to load. If you want to have the flexibibility to modify and run the model without loading the graphics, you may want to execute the plotting script separately.

### Creating a Model Script

Create a new file called `run_model.jl` in your project directory with the following content:

``` {.julia file=examples/model/alcap/run.jl}

```

The output will be saved in the `HDF5` format. This is a binary file and should not be tracked by version control. It can be accessed using CarboKitten's built in tools (see below for plotting) or standard libraries for reading HDF5 for Python or R. 

Please note that this script puts all new functions and variables in a module (`Module Script`). This is not obligatory, but recommended. You can [read more about Julia modules](https://docs.julialang.org/en/v1/manual/modules/) if you want to understand why this is useful.

### Running the Script

From your project directory, run:

```bash
julia --project=. run_model.jl
```

Alternatively, from the Julia REPL with the project activated:

```julia
include("run_model.jl")
```

The model will run and save output to `data/output/alcap-example.h5`. This may take a few minutes depending on your system.

## Visualizing Results

### Installing Visualization Dependencies

To visualize the model output, you need to install a Makie backend. For interactive plots, use GLMakie:

```julia
(MyCarboKittenProject) pkg> add GLMakie
```

This step will likely take longer than installing CarboKitten itself. On a slow computer, it can take more than 10 min. You will only have to do it once per project, although 

### Creating a Plotting Script

Create a new file called `plot_results.jl`:

```{.julia file=examples/model/alcap/plot_alcap.jl}
using GLMakie
using CarboKitten.Visualization

# Activate the GLMakie backend
GLMakie.activate!()

# Generate a summary plot
fig = summary_plot("my-first-model.h5")

# Display the plot
display(fig)

# Optionally save to file
save("my-first-model.png", fig)
```

### Running the Plotting Script

```bash
julia --project=. plot_results.jl
```

Or from the Julia REPL:

```julia
include("plot_results.jl")
```

The first time you run this, it may take a while to compile. Subsequent runs will be much faster.

## Understanding the Output

The summary plot shows:

- **Sediment profile**: Cross-section showing facies distribution and unconformities
- **Wheeler diagram**: Time-stratigraphic representation showing sedimentation rates and dominant facies
- **3D topography**: Three-dimensional view of the final platform geometry
- **Sea level curve**: The sea level variation through time
- **Production curves**: Carbonate production rates for each facies across water depths

## Next Steps

- Follow the [tutorial](first_tutorial.md)
- Explore different [input parameters](model-alcap.md) to customize your model
- Learn about [data export](data-export.md) options for further analysis
- Read about the [ALCAP model](model-alcap.md) components and theory
- Check out the [visualization](visualization.md) options for different plot types

## Multi-threading

To speed up model runs, use multiple threads:

```bash
julia --project=. -t 4 run_model.jl
```

This runs Julia with 4 threads. Adjust the number based on your CPU.

## Troubleshooting

If you encounter issues:

1. Make sure you're using Julia 1.10 or later: `julia --version`
2. Verify your project environment is activated (the REPL prompt should show your project name)
3. Try updating all packages: `pkg> up` in package mode
4. Check the [documentation](https://mindthegap-erc.github.io/CarboKitten.jl) for more details
5. Report bugs at [github.com/MindTheGap-ERC/CarboKitten.jl/issues](https://github.com/MindTheGap-ERC/CarboKitten.jl/issues) 
6. Ask for help at [https://github.com/MindTheGap-ERC/CarboKitten.jl/discussions](https://github.com/MindTheGap-ERC/CarboKitten.jl/discussions)
