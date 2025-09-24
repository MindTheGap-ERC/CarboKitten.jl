# Abstract Output

## Interface

``` {.julia file=src/Output/Abstract.jl}
module Abstract

export Data, DataColumn, DataSlice, DataVolume, Slice2, Header, DataHeader, Axes, AbstractOutput, Frame
export parse_multi_slice, data_kind, new_output, add_data_set, set_attribute, state_writer, frame_writer

using Unitful

const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}
const Amount = typeof(1.0u"m")
const Sediment = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")

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

## Run model

On top of this we have defined a `run_model` method that writes output in some form.

``` {.julia #run-model-output}
"""
    run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where M

Run a model and save the output to `output`.
"""
function run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where M
    state = M.initial_state(input)
    write_state = state_writer(input, output)
    write_frame = frame_writer(input, output)

    # create a group for every output item
    for (k, v) in input.output
        add_data_set(output, k, v)
    end
    write_state(1, state)

    run_model(Model{M}, input, state) do w, df
        # write_frame chooses to advance in a dataset
        # or just to increment on the current frame
        write_frame(w, df)
        # write_state only writes one in every write_interval
        # and does no accumulation
        write_state(w+1, state)
    end

    return output
end
```


``` {.julia file=src/Output/RunModel.jl}
module RunModel

import ...CarboKitten: run_model, Model
using ..Abstract 

<<run-model-output>>

end
```
