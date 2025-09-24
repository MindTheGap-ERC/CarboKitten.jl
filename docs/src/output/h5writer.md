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
function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where {M}
    output = H5Output(input, filename)
    run_model(Model{M}, input, output)
    return filename
end
```

## Implementation

``` {.julia file=src/Output/H5Writer.jl}
module H5Writer

using HDF5
using Unitful

import ...CarboKitten: run_model, Model, AbstractOutput, AbstractInput, OutputSpec, AbstractState

using ...CarboKitten: time_axis, box_axes
using ...Components.WaterDepth: initial_topography

using ...Utility: in_units_of
using ..Abstract
import ..Abstract: add_data_set, set_attribute, frame_writer, state_writer

mutable struct H5Output <: AbstractOutput
    header::Header
    fid::HDF5.File
end

function H5Output(input::AbstractInput, filename::String)
    t_axis = time_axis(input.time)
    x_axis, y_axis = box_axes(input.box)
    axes = Axes(x=x_axis, y=y_axis, t=t_axis)
    h0 = initial_topography(input)
    sl = input.sea_level.(t_axis)

    header = Header(
        tag=input.tag,
        axes=axes,
        Δt=input.time.Δt,
        time_steps=input.time.steps,
        grid_size=input.box.grid_size,
        n_facies=length(input.facies),
        initial_topography=h0,
        sea_level=sl,
        subsidence_rate=input.subsidence_rate,
        data_sets=Dict(),
        attributes=Dict())

    fid = h5open(filename, "w")
    create_group(fid, "input")

    finalizer(H5Output(header, fid)) do x
        close(x.fid)
    end
end

axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)

slice_str(::Colon) = ":"
slice_str(r::AbstractRange{Int}) = "$(r.start):$(r.stop)"
slice_str(i::Int) = "$(i)"

is_column(::Int, ::Int) = true
is_column(_, _) = false

function add_data_set(out::H5Output, name::Symbol, spec::OutputSpec)
    header = out.header
    fid = out.fid

    nf = header.n_facies
    nw = div(header.time_steps, spec.write_interval)
    size = axis_size.(spec.slice, header.grid_size)

    grp = HDF5.create_group(fid, string(name))
    attrs = attributes(grp)
    attrs["slice"] = join(slice_str.(spec.slice), ",")
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

function get_group(fid::HDF5.File, name::String)
    path = split(name, "/")
    gid = fid["input"]
    for n in path[1:end-1]
        if !haskey(gid, n)
            create_group(gid, n)
        end
        gid = gid[n]
    end
    return gid
end

function set_attribute(out::H5Output, name::String, value::AbstractArray{T,Dim}) where {T,Dim}
    gid = get_group(out.fid, name)
    gid[name] = value
end

function set_attribute(out::H5Output, name::String, value)
    gid = get_group(out.fid, name)
    attr = attributes(gid)
    attr[name] = value
end

function state_writer(input::AbstractInput, out::H5Output)
    output_spec = input.output
    fid = out.fid
    grid_size = out.header.grid_size

    function (idx::Int, state::AbstractState)
        for (k, v) in output_spec
            size = axis_size.(v.slice, grid_size)
            if mod(idx - 1, v.write_interval) == 0
                fid[string(k)]["sediment_thickness"][:, :, div(idx - 1, v.write_interval)+1] =
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
            tgt[:, :, :, div(idx - 1, v.write_interval)+1] +=
                (reshape(src[:, v.slice...], (n_f, size...)) |> in_units_of(u"m"))
        end

        for (k, v) in input.output
            n_writes = div(input.time.steps, v.write_interval)
            if div(idx-1, v.write_interval) + 1 <= n_writes
                grp = out.fid[string(k)]
                try_write(grp["production"], frame.production, v)
                try_write(grp["disintegration"], frame.disintegration, v)
                try_write(grp["deposition"], frame.deposition, v)
            end
        end
    end
end

<<hdf5-run-model>>

end
```

## Tests

```{.julia file=test/Output/H5WriterSpec.jl}
module H5WriterSpec

using CarboKitten
using HDF5
using CarboKitten.Output: write_frame, create_dataset, Frame
using CarboKitten.Output.H5Writer: H5Output
using Unitful
using Test

const DummyFacies = [
    ALCAP.Facies(
        viability_range = (0, 0),
        activation_range = (0, 0),
        maximum_growth_rate=0.0u"m/Myr",
        extinction_coefficient=0.0u"m^-1",
        saturation_intensity=0.0u"W/m^2",
        diffusion_coefficient=0.0u"m/yr")]

const input = ALCAP.Input(
    tag="test",
    box=Box{Periodic{2}}(grid_size=(5, 1), phys_scale=5.0u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=10),
    output=Dict(
        :wi1 => OutputSpec(slice=(:,:), write_interval=1),
        :wi2 => OutputSpec(slice=(:,:), write_interval=2),
        :wi3 => OutputSpec(slice=(:,:), write_interval=3),
        :wi4 => OutputSpec(slice=(:,:), write_interval=4)),
    ca_interval=1,
    initial_topography=(x, y) -> -0.0u"m",
    sea_level = t -> 0.0u"m",
    subsidence_rate=0.0u"m/Myr",
    disintegration_rate=0.0u"m/Myr",
    insolation=0.0u"W/m^2",
    sediment_buffer_size=0,
    depositional_resolution=0.0u"m",
    facies=DummyFacies)

const filename = "testH5.h5"


@testset "Components/H5Writer" begin

    mktempdir() do path
        fpath = joinpath(path, filename)
        h5open(fpath, "w") do fid
            create_group(fid, "input")
            ALCAP.write_header(fid, input)

            for (k, v) in input.output
                create_ck_group(fid, input, k, v)
            end

            # create a frame of ones to be the deposition etc. each time step
            inc = zeros(Frame, input)
            inc.production .= 1.0u"m"
            inc.deposition .= 1.0u"m"
            inc.disintegration .= 1.0u"m"

            for t = 1:input.time.steps
                write_frame(fid, input, t, inc)
            end
        end

        @testset "size of output array" begin
            h5open(fpath, "r") do f
                @test size(f["wi1"]["deposition"][])[4] == 10
                @test size(f["wi2"]["deposition"][])[4] == 5
                @test size(f["wi3"]["deposition"][])[4] == 3
                @test size(f["wi4"]["deposition"][])[4] == 2
            end
        end

        @testset "frame written only every write_interval" begin
            h5open(fpath, "r") do f
                @test all(f["wi1"]["deposition"][] .≈ attrs(f["wi1"])["write_interval"])
                @test all(f["wi2"]["deposition"][] .≈ attrs(f["wi2"])["write_interval"])
                @test all(f["wi3"]["deposition"][] .≈ attrs(f["wi3"])["write_interval"])
                @test all(f["wi4"]["deposition"][] .≈ attrs(f["wi4"])["write_interval"])
            end
        end
    end
end

end
```
