# HDF5 Output

We write output to HDF5.

``` {.julia file=src/Components/H5Writer.jl}
@compose module H5Writer
    using ..Common
    @mixin Boxes, TimeIntegration

    @kwdef struct DataFrame
        step::Int
        disintegration::Array{Amount,3}    # facies, x, y
        production::Array{Amount,3}
        deposition::Array{Amount,3}
        sediment_height::Array{Amount,2}
    end

    function create_dataset(fid, input::AbstractInput)
        nf = n_facies(input)
        HDF5.create_dataset(fid, "production", datatype(Float64),
            dataspace(nf, input.box.grid_size..., input.time.steps),
            chunk=(nf, input.box.grid_size..., 1))
        HDF5.create_dataset(fid, "disintegration", datatype(Float64),
            dataspace(nf, input.box.grid_size..., input.time.steps),
            chunk=(nf, input.box.grid_size..., 1))
        HDF5.create_dataset(fid, "deposition", datatype(Float64),
            dataspace(nf, input.box.grid_size..., input.time.steps),
            chunk=(nf, input.box.grid_size..., 1))
        HDF5.create_dataset(fid, "sediment_height", datatype(Float64),
            dataspace(input.box.grid_size..., input.time.steps),
            chunk=(input.box.grid_size..., 1))
    end

    function write_frame(fid, frame::DataFrame)
        step = frame.step
        fid["production"][:, :, :, step] = frame.production |> in_units_of(u"m")
        fid["disintegration"][:, :, :, step] = frame.disintegration |> in_units_of(u"m")
        fid["deposition"][:, :, :, step] = frame.deposition |> in_units_of(u"m")
        fid["sediment_height"][:, :, step] = frame.sediment_height |> in_units_of(u"m")
    end

    function run(::Type{Model}, input::AbstractInput)
    end
end
```
