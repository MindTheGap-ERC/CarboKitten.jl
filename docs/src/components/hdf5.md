# HDF5 Output

```component-dag
CarboKitten.Components.H5Writer
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

``` {.julia #hdf5-output-spec}
@kwdef struct Input <: AbstractInput
    output = Dict(:full => OutputSpec((:,:), 1)) 
end
```

The default is to write all output, which is fine for smaller runs. The `slice` argument of `OutputSpec` will accept three different forms:

- `(:, :)` (default) output the full area of the model.
- `(<n>, :)` or `(:, <n>)` output a slice of the model, either with a fixed $x$ or a fixed $y$ coordinate. In our examples we always have the $x$ axis orthogonal to the shoreline, so slicing with a fixed $y$ (the second form) is what we use.
- `(<m>, <n>)`, output a column of the model. If you have a very precise experiment workflow, this could be of use. You'll have to specify each column as a separate output. Most of the time though, we can extract columns from slice data in a post-processing stage, so if all your columns have the same $y$ coordinate, taking a slice is the preferred option.

## Main loop
The `H5Writer` component runs through an overloaded version of `run_model`. For example:

```julia
run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
```

``` {.julia #hdf5-run-model}
"""
    run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M

Run a model and write output to HDF5. Here `M` should be a model, i.e. a
module with `initial_state`, `step!` and `write_header` defined. Example:

    run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
"""
function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M        
    state = M.initial_state(input)

    h5open(filename, "w") do fid
        create_group(fid, "input")
        M.write_header(fid, input)

        # create a group for every output item
        for (k, v) in input.output
            create_ck_group(fid, input, k, v)
        end
        write_state(fid, input, 1, state)

        run_model(Model{M}, input, state) do w, df
            # write_frame chooses to advance in a dataset
            # or just to increment on the current frame
            write_frame(fid, input, w, df)
            # write_state only writes one in every write_interval
            # and does no accumulation
            write_state(fid, input, w+1, state)
        end
    end

    return filename
end
```

## Implementation

``` {.julia file=src/Components/H5Writer.jl}
@compose module H5Writer
    using ..Common
    import ..Models: Frame

    using HDF5

    import ...CarboKitten: run_model, Model
    import ...OutputData: AbstractOutputSpec

    @mixin Boxes, TimeIntegration, FaciesBase, WaterDepth

    <<hdf5-output-spec>>

	axis_size(::Colon, a::Int) = a
	axis_size(::Int, _) = 1
	axis_size(r::AbstractRange{Int}, _) = length(r)

	slice_str(::Colon) = ":"
	slice_str(r::AbstractRange{Int}) = "$(r.start):$(r.stop)"
	slice_str(i::Int) = "$(i)"

    is_column(::Int, ::Int) = true
    is_column(_, _) = false

    function create_ck_group(fid, input::AbstractInput, name::Symbol, spec::AbstractOutputSpec)
        nf = n_facies(input)
        nw = div(input.time.steps, spec.write_interval)

        size = axis_size.(spec.slice, input.box.grid_size)

        grp = HDF5.create_group(fid, string(name))
        attrs = attributes(grp)
		attrs["slice"] = join(slice_str.(spec.slice),",")
        attrs["write_interval"] = spec.write_interval

        HDF5.create_dataset(grp, "production", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)	
        HDF5.create_dataset(grp, "disintegration", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)
        HDF5.create_dataset(grp, "deposition", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)
        HDF5.create_dataset(grp, "sediment_thickness", datatype(Float64),
            dataspace(size..., nw + 1),
            chunk=(size..., 1), deflate=3)
    end

    function write_state(fid, input::AbstractInput, idx::Int, state::AbstractState)
        for (k, v) in input.output
			size = axis_size.(v.slice, input.box.grid_size)
            if mod(idx-1, v.write_interval) == 0
                fid[string(k)]["sediment_thickness"][:, :, div(idx-1, v.write_interval)+1] = 
                    (is_column(v.slice...) ?
                        Float64[state.sediment_height[v.slice...] |> in_units_of(u"m");] :
                        reshape(state.sediment_height[v.slice...], size) |> in_units_of(u"m"))
            end
        end
    end

    function write_frame(fid, input::AbstractInput, idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, v) = ()
		n_f = n_facies(input)
        function try_write(tgt, src::AbstractArray, v)
			size = axis_size.(v.slice, input.box.grid_size)
            tgt[:, :, :, div(idx-1, v.write_interval) + 1] += 
				(reshape(src[:, v.slice...], (n_f, size...)) |> in_units_of(u"m"))
        end

        for (k, v) in input.output
            n_writes = div(input.time.steps, v.write_interval)
            # only write if still within bounds of array
            if div(idx-1, v.write_interval)+1 <= n_writes
                grp = fid[string(k)]
                try_write(grp["production"], frame.production, v)
                try_write(grp["disintegration"], frame.disintegration, v)
                try_write(grp["deposition"], frame.deposition, v)
            end
        end
    end

    <<hdf5-run-model>>
end
```
