# Output

## Data Structures

``` {.julia file=src/OutputData.jl}
module OutputData

export Data, DataColumn, DataSlice, DataVolume, Slice2, Header, DataHeader, Axes, AbstractOutput, Frame
export parse_multi_slice, data_kind, new_output, add_data_set, set_attribute, state_writer, frame_writer

using Unitful

const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}
const Amount = typeof(1.0u"m")
const Sediment = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")

abstract type AbstractOutputSpec end
abstract type AbstractInput end
abstract type AbstractOutput end
abstract type AbstractState end

@kwdef struct OutputSpec <: AbstractOutputSpec
    slice::Slice2 = (:, :)
    write_interval::Int = 1
end

@kwdef struct Axes
    x::Vector{Length}
    y::Vector{Length}
    t::Vector{Time}
end

@kwdef struct DataHeader
    kind::Symbol
    slice::Slice2
    write_interval::Int
end

@kwdef struct Header
    tag::String
    axes::Axes

    Δt::Time
    time_steps::Int
    grid_size::NTuple{2,Int}
    n_facies::Int

    initial_topography::Matrix{Amount}
    sea_level::Vector{Length}
    subsidence_rate::Rate
    data_sets::Dict{Symbol, DataHeader}
    attributes::Dict{Symbol, Any} = Dict()
end

@kwdef struct Data{F, D}
	slice::Slice2
	write_interval::Int
	# Julia doesn't allow to say Array{Amount,D+1} here
    disintegration::Array{Amount,F}
	production::Array{Amount,F}
    deposition::Array{Amount,F}
    sediment_thickness::Array{Amount,D}
end

const DataVolume = Data{4, 3}
const DataSlice = Data{3, 2}
const DataColumn = Data{2, 1}

@kwdef struct Frame
    disintegration::Union{Array{Sediment,3},Nothing} = nothing   # facies, x, y
    production::Union{Array{Sediment,3},Nothing} = nothing
    deposition::Union{Array{Sediment,3},Nothing} = nothing
end

count_ints(::Int, args...) = 1 + count_ints(args...)
count_ints(_, args...) = count_ints(args...)
count_ints() = 0

reduce_slice(s::Tuple{Colon, Colon}, x, y) = (x, y)
reduce_slice(s::Tuple{Int, Colon}, y::Int) = (s[1], y)
reduce_slice(s::Tuple{Colon, Int}, x::Int) = (x, s[2])

Base.getindex(v::Data{F,D}, args...) where {F, D} = let k = count_ints(args...)
	Data{F-k, D-k}(
		reduce_slice(v.slice, args...),
		v.write_interval,
		v.disintegration[:, args..., :],
		v.production[:, args..., :],
		v.deposition[:, args..., :],
		v.sediment_thickness[args..., :])
end

function parse_slice(s::AbstractString)
	if s == ":"
		return (:)
	end

	elements =  split(s, ":")
	if length(elements) == 1
		return parse(Int, s)
	end

	a, b = elements
	return parse(Int, a):parse(Int, b)
end

parse_multi_slice(s::AbstractString) = Slice2(parse_slice.(split(s, ",")))

data_kind(::Int, ::Int) = :column
data_kind(::Int, _) = :slice
data_kind(_, ::Int) = :slice
data_kind(_, _) = :volume
data_kind(spec::OutputSpec) = data_kind(spec.slice...)

"""
    new_output(::Type{T}, input)

Create a new output object of type `T`, given `input`.
"""
function new_output end

"""
    add_data_set(out::T, name::Symbol, spec::OutputSpec)

Add a data set to the output object.
"""
function add_data_set end

"""
    set_attribute(out::T, name::Symbol, value::Any)

Set an attribute in the output object.
"""
function set_attribute end

"""
    write_sediment_thickness(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

Write the sediment thickness to the output object. The `idx` should be corrected for write
interval. That is, `idx` should range from `1` to `n_writes` for the named data set. This
function should be implemented for 0, 1, and 2 dimensional arrays, corresponding to writing
column, slice or volume data.

If your output object type doesn't conform to the standard CarboKitten data layout, you may
choose to not implement this function and implement `state_writer` and `frame_writer` instead.
The same goes for `write_production`, `write_disintegration` and `write_deposition`.
"""
function write_sediment_thickness end

"""
    write_production(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_production end

"""
    write_disintegration(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_disintegration end

"""
    write_deposition(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_deposition end

"""
    state_writer(input::AbstractInput, out::T)

Returns a `function (idx::Int, state::AbstractState)`.

Write the state of the simulation to the output object. It is the responsibility
of the implementation to choose to write or not based on the set write interval.

The default implementation writes the state for all output data sets, and calls
`write_sediment_thickness`.
"""
function state_writer(input::AbstractInput, out)
    output_sets = input.output
    grid_size = input.box.grid_size

    return function(idx::Int, state::AbstractState)
        for (k, v) in output_sets
            if mod(idx-1, v.write_interval) == 0
                write_sediment_thickness(
                    out, k, div(idx-1, v.write_interval)+1,
                    view(state.sediment_height, v.slice...)) 
            end
        end
    end
end

"""
    frame_writer(input::AbstractInput, out::T)

Returns a `function (idx::Int, state::AbstractState)`.

Write the state of the simulation to the output object. It is the responsibility
of the implementation to choose to write or not based on the set write interval.

The default implementation writes the state for all output data sets, and calls
`write_sediment_thickness`.
"""
function frame_writer(input::AbstractInput, out)
    n_f = length(input.facies)
    grid_size = input.box.grid_size

    return function(idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, k, v) = ()
        function try_write(write::F, src, k, v) where {F}
            write(out, k, div(idx-1, v.write_interval) + 1,
                  view(src, :, v.slice...))
        end

        for (k, v) in input.output
            try_write(write_production, frame.production, k, v)
            try_write(write_disintegration, frame.disintegration, k, v)
            try_write(write_deposition, frame.deposition, k, v)
        end
    end
end

end
```

## Memory Writer

Sometimes we don't want to write to HDF5, but just get a `DataVolume` directly.

``` {.julia file=src/MemoryWriter.jl}
module MemoryWriter

using ..OutputData
import ..OutputData:
    new_output, add_data_set, set_attribute, write_sediment_thickness,
    write_production, write_disintegration, write_deposition, AbstractOutputSpec

using ..Components.Common
using ..Components.WaterDepth: initial_topography
using ..CarboKitten: time_axis, box_axes

struct MemoryOutput <: AbstractOutput
    header::Header
    data_volumes::Dict{Symbol, DataVolume}
    data_slices::Dict{Symbol, DataSlice}
    data_columns::Dict{Symbol, DataColumn}
end

function new_output(::Type{MemoryOutput}, input::AbstractInput)
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

    return MemoryOutput(header, Dict(), Dict(), Dict())
end

axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)

function add_data_set(out::MemoryOutput, label::Symbol, spec::AbstractOutputSpec)
    h = DataHeader(data_kind(spec), spec.slice, spec.write_interval)
    out.header.data_sets[label] = h

    slice = spec.slice
    write_interval = spec.write_interval

    full_size = out.header.grid_size
    n_steps = div(out.header.time_steps, write_interval)
    n_facies = out.header.n_facies

    if h.kind == :volume
        size = axis_size.(slice, full_size)
        out.data_volumes[label] = DataVolume(
            slice, write_interval,
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, size..., n_steps + 1))
    elseif h.kind == :slice
        size = axis_size.(slice, full_size)
        slice_size = size[1] == 1 ? size[2] : size[1]
        out.data_slices[label] = DataSlice(
            slice, write_interval,
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, slice_size, n_steps + 1)) 
    elseif h.kind == :column
        out.data_columns[label] = DataColumn(
            slice, write_interval,
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_steps + 1))
    end
end

function set_attribute(out::MemoryOutput, name::Symbol, value::Any)
    out.header.attributes[name] = value
end

write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 0}) =
    out.data_columns[label].sediment_thickness[idx] .= data
write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_slices[label].sediment_thickness[:, idx] .= data
write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_volumes[label].sediment_thickness[:, :, idx] .= data

write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].production[:, idx] .+= data
write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].production[:, :, idx] .+= data
write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].production[:, :, :, idx] .+= data

write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].disintegration[:, idx] .+= data
write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].disintegration[:, :, idx] .+= data
write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].disintegration[:, :, :, idx] .+= data

write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].deposition[:, idx] .+= data
write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].deposition[:, :, idx] .+= data
write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].deposition[:, :, :, idx] .+= data

end
```
