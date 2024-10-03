# Introduction to CarboKitten.jl

## Motivation and Aim
In completeness of stratigraphy would substantially bias our interpretations on the archives. However, if ther

This project is written with `Julia`, an open-source and free language. For more information, please click [here](https://docs.julialang.org/en/v1/)

## Install Julia
Details of downloading and installing `Julia` could be find [here](https://julialang.org/downloads/). For example, if you are using Windows, please type the following command in command prompt:

```winget install julia -s msstore``` 

Alternatively, the you could use `juliaup` from [here](github.com/JuliaLang/juliaup) to download julia.

## REPL (read-eval-print loop)
REPL allows users to interactively use Julia. That is to say, REPL takes your inputs, executes them, and returns the result to the you. 

### Easiest way
You can simply click julia.exe to start a REPL.

### VSCode
VScode is an integrated development environment that supports many mainstream language. You can start a julia REPL in VSCode. 

Download VSCode from https://code.visualstudio.com/. after installing it, find 'extension' tab at the left bar, and type in julia, added this extension. The it should good to go. You can start a julia REPL by `Ctrl + Shft + P`. 

### Pluto (I think for workshop, this is the best option)
If you don't want to install too many things on your computer, you can try Pluto notebook, similar to Jupyter notebook.

enter package mode by:

 ```]```

```add Pluto```

```import Pluto```

```Pluto.run()```

The browser will pop up and you can create a new notebook.

Make sure you are in the path of cloned repo for Pluto.

## Clone the repo
Go to the folder you would like to work on.

```cd ./your/path```

You can clone the repository from github in powershell:

```git clone https://github.com/MindTheGap-ERC/CarboKitten.jl.git ```

Next step is to activate environment in Julia REPL (in Package mode):

```activate .``` or ```activate ../workenv```

Install all dependencies (still in Package mode):

```instantiate```

Now you are ready to go.

## Run model Example
There are few examples were ready for you. 

```include("./examples/production-only/ca-slope.jl")```

You can type in julia REPL:

```readdir("./examples")```

to see how many examples we have. You can try to run all examples by changing the `$name.jl`. The details of the examples could be found in xxx.md. 

The results are saved in "./data" as HDF5 file. (we have to correspond the example name before workshop)

## Visualisation and data extraction
Visualisation are not included in CarboKitten.jl, and are implemented in a Julia package extension.

Better to have a ready-to-use Pluto notebook?

- Visualisation of cross-section and wheeler diagram
Need `Makie`.

in package model, `add GLMakie`.

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

- export data
You can copy paster the following code in Pluto to export data at certain point. In this case, it exported 10, 30, 50 and 70 from the result.

```
using CarboKitten.Export: CSV, data_export

data_export(
    CSV(tuple.(10:20:70, 25),
      :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
      :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
      :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
      :metadata => "$(PATH)/$(TAG).toml"),
    "$(PATH)/$example.h5")
```

- Stratigraphic columns.

- Age-depth model 

## Change input
You could copy one example in Pluto, and try to change the input values and run the new example.

We will take sea-level as an example. The default input sea-level curve -> 0, and it means that sea-level did not change.

You can define your own sea-level curve. For example, replace `sea_level = 0.0u"m"` to `sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD)u"m"`. Where `AMPLITUDE` is the amplitude of sea-level and `PERIOD` is period of the sea-level cycle.

You can also replace the sea-level curve with real data. Herein, we prepare curve from [miller_phanerozoic_2005](@cite) in "./data/millers.csv". This is the observed sea-level curve from last 2Ma, with timestep of 1kyr. You can type this in Pluto:

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

`sea_level = t -> sealevel_curve(t,filepath)`

Remember to add `DataFrames`, `CSV`, `Interpolations` Packages, and change the output name accordingly. 