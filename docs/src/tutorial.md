# Introduction to CarboKitten.jl

## Motivation and Aim
In completeness of stratigraphy would substantially bias our interpretations on the archives. However, if ther

This project is written with `Julia`, an open-source and free language. For more information, please click [here](https://docs.julialang.org/en/v1/)

## Install Julia
Details of downloading and installing `Julia` could be find [here](https://julialang.org/downloads/). For example, if you are using Windows, please type the following command in command prompt:

```shell
winget install julia -s msstore
``` 

Alternatively, the you could use `juliaup` from [here](github.com/JuliaLang/juliaup) to download julia.

## REPL (read-eval-print loop)
REPL allows users to interactively use Julia. That is to say, REPL takes your inputs, executes them, and returns the result to the you. 

### Easiest way
You can simply click julia.exe to start a REPL.

### VSCode
VScode is an integrated development environment that supports many mainstream language. You can start a julia REPL in VSCode. 

Download VSCode from https://code.visualstudio.com/. after installing it, find 'extension' tab at the left bar, and type in julia, added this extension. The it should good to go. You can start a julia REPL by `Ctrl + Shft + P`. 

### Pluto (We will stick to Pluto in the workshop)
If you don't want to install too many things on your computer, you can try Pluto notebook, similar to Jupyter notebook.

enter package mode by:

```juliarepl
julia> using Pluto
# answer question [y]
julia> Pluto.run()
```

The browser will pop up and you can create a new notebook. To stop Pluto, please enter `Ctrl + C`.

You can run code blocks by clicking the small triangle at the right down side of each block.

## Clone the repo
Go to the folder you would like to clone and store the repo locally (in powershell).

```cd ./your/path```

You can clone the repository from github in powershell:

```git clone https://github.com/MindTheGap-ERC/CarboKitten.jl.git ```

```cd CarboKitten.jl```

Here you can enter `ls` or `tree.` to see the structure of this repo.

Next step is to open julia REPL (as described in above session), and direct to the folder you are working on:

```cd("./your/working/folder")```

Then we have to activate environment in Julia REPL (in Package mode)

```activate .``` 

The `.` here means your project is activated in current directory.

Install all dependencies (still in Package mode):

```instantiate```

Now you are ready to go.

(I think this step is better for people to use if we navigate them to the existing notebook in "./notebooks/examplenotebook.jl"), where all the code is available.

## Run model Example
There are few examples were ready for you. 

```include("../examples/production-only/ca-slope.jl")```

You can type in julia REPL:

```readdir("../examples")```

to see how many examples we have. You can try to run all examples by changing the `$name.jl`. The details of the examples could be found in xxx.md. 

The results are saved in "./data" as HDF5 file. (we have to correspond the example name before workshop)

## Visualisation and data extraction

(Better to have a ready-to-use Pluto notebook?

### Visualisation of cross-section and wheeler diagram

`using GLMakie`.

Visualise wheeler diagram:

```using CarboKitten.Visualization: sediment_profile, wheeler_diagram```

```
module Script

function main()
  header, data = read_slice("data/$example.h5", :, 25)
  fig = wheeler_diagram(header, data)
  save("docs/src/_fig/e$xample.png", fig)
end

end

Script.main()
```

Visualise cross-section:

```
  fig = Visualization.sediment_profile("data/$example.h5", 25)
  save("docs/src/_fig/$example.png", fig)
```

### export data
You can copy paster the following code in Pluto to export data at certain point. Herein, in `CSV(tuple.(50,25))`, the `50` means it's 50 km away from the shoal. You can vary this value between 1 and 100 and observe the changes.

```
using CarboKitten.Export: CSV, data_export

data_export(
    CSV(tuple.(50, 25),
      :sediment_accumulation_curve => "$example_sac.csv",
      :age_depth_model => "$example_adm.csv",
      :stratigraphic_column => "$example_sc.csv",
      :metadata => "$example.toml"),
    "$example.h5")
```

The sediment-accumulation-curve, age-depth model and stratigraphy-columns should now be saved as csv files.

### Stratigraphic columns

```
import CarboKitten.Visualization: stratigraphic_column!
using CarboKitten.Export: stratigraphic_column, age_depth_model
```

```
wait for the function is ready
```

### Age-depth model 

```
wait for the function is ready
```

## Change input
You could copy one example into Pluto, and try to change the input values and run the new example.

We will take changing sea-level for a instance. The default input sea-level curve -> 0, and it means that sea-level did not change.

You can define your own sea-level curve. For example, replace `sea_level = 0.0u"m"` to `sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD)u"m"`. Where `AMPLITUDE` is the amplitude of sea-level and `PERIOD` is period of the sea-level cycle.

You can also replace the sea-level curve with observed data. Herein, we prepare curve from [miller_phanerozoic_2005](@cite) in "./data/millers.csv". This is the observed sea-level curve from last 2Ma, with timestep of 1kyr. You can type this in Pluto:
Please first add `DataFrames`, `CSV`, `Interpolations` Packages (in package mode, type `add CSV`)

```
using DataFrames
using CSV
using Interpolations
function sealevel_curve(t,filepath)
    data = DataFrame(CSV.File(filepath))
    data = hcat(collect(0:length(data.sl).-1)./1000, data.sl) #The output sealevel curve from R does not have time tab and have to add it 
    x = linear_interpolation(data[:,1], data[:,2])
    return x(t)
end 
```

`sea_level = t -> sealevel_curve(t,filepath) .* 1.0u"m"`

Please remember to change the output HDF5 name accordingly. 