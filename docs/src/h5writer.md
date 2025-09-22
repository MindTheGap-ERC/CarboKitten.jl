# HDF5 Writer

## Main loop
The `H5Writer` runs through an overloaded version of `run_model`. For example:

```julia
run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
```

``` {.julia #hdf5-run-model}
"""
    run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M

Run a model and write output to HDF5.

    run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
"""
function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M
    output = H5Output(input, filename)
    run_model(Model{M}, input, output)
    return filename
end
```

## Implementation

``` {.julia file=src/Output/H5Writer.jl}
module H5Writer

using HDF5

import ...CarboKitten: run_model, Model
import ..OutputData: Frame, AbstractInput, AbstractOutput, OutputSpec, AbstractState,
                     add_data_set, set_attribute, frame_writer, state_writer

struct H5Output <: AbstractOutput
    header::Header
    fid::HDF5.File
end

function H5Output(input::AbstractInput, filename::String)
    t_axis = time_axis(input.time)
    x_axis, y_axis = box_axes(input.box)
    axes = Axes(x = x_axis, y = y_axis, t = t_axis)
    h0 = initial_topography(input)
    sl = input.sea_level.(t_axis)

    header = Header(
        tag = input.tag,
        axes = axes,
        Δt = input.time.Δt,
        time_steps = input.time.steps,
        grid_size = input.box.grid_size,
        n_facies = length(input.facies),
        initial_topography = h0,
        sea_level = sl,
        subsidence_rate = input.subsidence_rate,
        data_sets = Dict(),
        attributes = Dict())

    fid = HDF5.open(filename)
    return H5Output(header, fid)
end

axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)

slice_str(::Colon) = ":"
slice_str(r::AbstractRange{Int}) = "$(r.start):$(r.stop)"
slice_str(i::Int) = "$(i)"

is_column(::Int, ::Int) = true
is_column(_, _) = false

function add_data_set(h5out::H5Output, name::Symbol, spec::OutputSpec)
    header = h5out.header
    fid = h5out.fid

    nf = header.n_facies
    nw = div(header.time_steps, spec.write_interval)
    size = axis_size.(spec.slice, header.grid_size)

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

function set_attribute(out::H5Output, name::Symbol, value::AbstractArray{T, Dim}) where { T, Dim }
    gid = fid["input"]
    gid[string(name)] = value
end

function set_attribute(out::H5Output, name::Symbol, value)
    gid = fid["input"]
    attr = attributes(gid)
    attr[string(name)] = value
end

function state_writer(input::AbstractInput, out::H5Output)
    output_spec = input.output
    fid = out.fid
    grid_size = out.header.grid_size

    function (idx::Int, state::AbstractState)
        for (k, v) in output_spec
            size = axis_size.(v.slice, grid_size)
            if mod(idx-1, v.write_interval) == 0
                fid[string(k)]["sediment_thickness"][:, :, div(idx-1, v.write_interval)+1] =
                    (is_column(v.slice...) ?
                        Float64[state.sediment_height[v.slice...] |> in_units_of(u"m");] :
                        reshape(state.sediment_height[v.slice...], size) |> in_units_of(u"m"))
            end
        end
    end
end

function frame_writer(input::AbstractInput, out::H5Output)
    n_f = out.header.n_facies
    grid_size = out.header.grid_size
    output_spec = input.output

    function (idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, v) = ()
        function try_write(tgt, src::AbstractArray, v)
            size = axis_size.(v.slice, grid_size)
            tgt[:, :, :, div(idx-1, v.write_interval) + 1] +=
                (reshape(src[:, v.slice...], (n_f, size...)) |> in_units_of(u"m"))
        end

        for (k, v) in input.output
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
