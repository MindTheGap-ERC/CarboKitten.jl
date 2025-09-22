# Output

```component-dag
CarboKitten.Components.Output
```

We write output to HDF5. In the `Input` struct the user can specify a dictionary of `OutputSpec`, specifying how much and at which interval to write output. Typically, you'd want a full topographic output with lower time resolution, and choose a transect with full time resolution. For example, on the ALCAP model:

```julia
const INPUT = ALCAP.Input(
    box = Box{Coast}(grid_size = (300, 150), phys_scale = 50.0u"m"),
    time = TimeProperties(Δt = 50.0u"yr", steps = 20_000),
    output = Dict(
        :topography => OutputSpec(write_interval = 200),
        :profile => OutputSpec(slice = (:, 75))),

    ...)  # add more options
```

Saving the full output of this simulation would take several hundreds of gigabytes, not gargantuan, but a bit unwieldy if you want to save many simulation runs. With this output specification, we cut down on this significantly.

``` {.julia #output-spec}
@kwdef struct Input <: AbstractInput
    output = Dict(:full => OutputSpec((:,:), 1))
end
```

The default is to write all output, which is fine for smaller runs. The `slice` argument of `OutputSpec` will accept three different forms:

- `(:, :)` (default) output the full area of the model.
- `(<n>, :)` or `(:, <n>)` output a slice of the model, either with a fixed $x$ or a fixed $y$ coordinate. In our examples we always have the $x$ axis orthogonal to the shoreline, so slicing with a fixed $y$ (the second form) is what we use.
- `(<m>, <n>)`, output a column of the model. If you have a very precise experiment workflow, this could be of use. You'll have to specify each column as a separate output. Most of the time though, we can extract columns from slice data in a post-processing stage, so if all your columns have the same $y$ coordinate, taking a slice is the preferred option.

``` {.julia file=src/Components/Output.jl}
@compose module Output
    <<output-spec>>
end
```
