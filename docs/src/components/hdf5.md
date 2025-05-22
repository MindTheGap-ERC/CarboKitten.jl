# HDF5 Output

```component-dag
CarboKitten.Components.H5Writer
```

We write output to HDF5.

``` {.julia file=src/Components/H5Writer.jl}
@compose module H5Writer
    using ..Common
    import ..Models: Frame

    using HDF5

    import ...CarboKitten: run_model, Model

    @mixin Boxes, TimeIntegration, FaciesBase, WaterDepth

    export run

    const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}

    @kwdef struct OutputSpec
        slice::Slice2 = (:, :)
        write_interval::Int = 1
    end

    @kwdef struct Input <: AbstractInput
        output = Dict(:full => OutputSpec((:,:), 1)) 
    end

    function create_ck_group(fid, input::AbstractInput, name::Symbol, spec::OutputSpec)
        nf = n_facies(input)
        nw = div(input.time.steps, spec.write_interval)

        axis_size(::Colon, a::Int) = a
        axis_size(::Int, _) = 1
        axis_size(r::AbstractRange{Int}, _) = length(r)

        size = axis_size.(spec.slice, input.box.grid_size)

        grp = HDF5.create_group(fid, string(name))
        attrs = attributes(grp)
        attrs["slice"] = string(spec.slice)
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
        HDF5.create_dataset(grp, "sediment_height", datatype(Float64),
            dataspace(size..., nw + 1),
            chunk=(size..., 1), deflate=3)
    end

    function write_state(fid, input::AbstractInput, idx::Int, state::AbstractState)
        for (k, v) in input.output
            if mod(idx, v.write_interval) == 0
                fid[string(k)]["sediment_height"][:, :, div(idx, v.write_interval)] = 
                    state.sediment_height[v.slice...] |> in_units_of(u"m")
            end
        end
    end

    function write_frame(fid, input::AbstractInput, idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, v) = ()
        try_write(tgt, src::AbstractArray, v) = if mod(idx, v.write_interval) == 0
            tgt[:, :, :, div(idx, v.write_interval)] = (src[:, v.slice...] |> in_units_of(u"m"))
        end

        for (k, v) in input.output
            grp = fid[string(k)]
            try_write(grp["production"], frame.production, v)
            try_write(grp["disintegration"], frame.disintegration, v)
            try_write(grp["deposition"], frame.deposition, v)
        end
    end

    """
        run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M

    Run a model and write output to HDF5. Here `M` should be a model, i.e. a
    module with `initial_state`, `step!` and `write_header` defined. Example:

    ```julia
    run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
    ```
    """
    function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M        
        state = M.initial_state(input)

        h5open(filename, "w") do fid
            create_group(fid, "input")
            M.write_header(fid, input)

            for (k, v) in input.output
                create_ck_group(fid, input, k, v)
            end
            write_state(fid, input, 1, state)

            run_model(Model{M}, input, state) do w, df 
                write_frame(fid, input, w, df)
                write_state(fid, input, w+1, state)
            end
        end

        return filename
    end
end
```
