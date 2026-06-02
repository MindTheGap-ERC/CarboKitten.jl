# ~/~ begin <<docs/src/output/abstract.md#src/Output/Abstract.jl>>[init]
module Abstract

import ...CarboKitten: set_attribute
import ...Algorithms: stratigraphic_column!

# Pull `cumulative_subsidence` from the Subsidence module so the methods we
# add here for `Header` extend the same generic function.
import CarboKitten.Components.Subsidence: cumulative_subsidence, deserialize_modifier, AbstractSubsidenceModifier

export Data, DataColumn, DataSlice, DataVolume, Slice2, Header, DataHeader, Axes, AbstractOutput, Frame
export parse_multi_slice, data_kind, new_output, add_data_set, set_attribute, state_writer, frame_writer
export cumulative_subsidence

using Unitful
using ...CarboKitten: OutputSpec, AbstractInput, AbstractState
using .Iterators: repeated

const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Slice2 = NTuple{2,Union{Int,Colon,UnitRange{Int}}}
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
    subsidence_rate_map::Union{Matrix{Rate},Nothing} = nothing #(optional, back-compatible): full per-cell rate map.
    #`nothing` means subsidence was uniform — use the scalar `subsidence_rate` instead.
    subsidence_modifiers::Vector{Any} = Any[] # modifier descriptions. Entries may be either live
        # `AbstractSubsidenceModifier` objects (when the header was built in memory)
        # or dict descriptors (when read from HDF5); `cumulative_subsidence` handles
        # both via `deserialize_modifier`.
    data_sets::Dict{Symbol,DataHeader}
    attributes::Dict{String,Any} = Dict()
end

@kwdef struct Data{F,D}
    slice::Slice2
    write_interval::Int
    # Julia doesn't allow to say Array{Amount,D+1} here
    disintegration::Array{Amount,F}
    production::Array{Amount,F}
    deposition::Array{Amount,F}
    sediment_thickness::Array{Amount,D}
    active_layer::Union{Array{Amount,F}, Nothing} = nothing
end

const DataVolume = Data{4,3}
const DataSlice = Data{3,2}
const DataColumn = Data{2,1}

@kwdef struct Frame
    disintegration::Union{Array{Sediment,3},Nothing} = nothing   # facies, x, y
    production::Union{Array{Sediment,3},Nothing} = nothing
    deposition::Union{Array{Sediment,3},Nothing} = nothing
end

count_ints(::Int, args...) = 1 + count_ints(args...)
count_ints(_, args...) = count_ints(args...)
count_ints() = 0

reduce_slice(s::Tuple{Colon,Colon}, x, y) = (x, y)
reduce_slice(s::Tuple{Int,Colon}, y::Int) = (s[1], y)
reduce_slice(s::Tuple{Colon,Int}, x::Int) = (x, s[2])

Base.getindex(v::Data{F,D}, args...) where {F,D} =
    let k = count_ints(args...)
        Data{F - k,D - k}(
            reduce_slice(v.slice, args...),
            v.write_interval,
            v.disintegration[:, args..., :],
            v.production[:, args..., :],
            v.deposition[:, args..., :],
            v.sediment_thickness[args..., :],
            v.active_layer == nothing ? nothing : v.active_layer[:, args..., :])
    end

function parse_slice(s::AbstractString)
    if s == ":"
        return (:)
    end

    elements = split(s, ":")
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
    stratigraphic_column(data)

Given a data set, compute the stratigrahpic column.
"""
function stratigraphic_column(data::Data{F, D}) where {F, D}
    net_deposition = data.deposition .- data.disintegration
    for c in eachslice(net_deposition, dims=(1:D...,))
        stratigraphic_column!(c)
    end
    return net_deposition
end

# =============================================================================
# Cumulative-subsidence helpers for post-hoc analysis
# =============================================================================

"""
    cumulative_subsidence(header::Header) -> (t -> Matrix{Length})

Closure form: returns a function that, given an absolute time `t`, returns the
per-cell cumulative subsidence (from `header.axes.t[1]`) as a `Matrix{Length}`
of shape `header.grid_size`.

Handles three cases:
  * Uniform rate, no modifiers — closed form `rate * (t - t0)`.
  * Per-cell map, no modifiers — `rate_map .* (t - t0)`.
  * With modifiers — piecewise integration via `CarboKitten.Subsidence`.

Modifiers in `header.subsidence_modifiers` may be either live modifier objects
or dict descriptors (as read from HDF5); both are accepted.
"""
function cumulative_subsidence(header::Header)
    t0 = header.axes.t[1]
    base = something(header.subsidence_rate_map,
                     fill(header.subsidence_rate, header.grid_size...))
    mods = AbstractSubsidenceModifier[deserialize_modifier(m) for m in header.subsidence_modifiers]
    # Delegate to the pure (base, modifiers, axes, t0) method
    return cumulative_subsidence(base, mods, header.axes.x, header.axes.y, t0)
end

"""
    cumulative_subsidence(header::Header, t::Time) -> Matrix{Length}

Direct evaluation: per-cell cumulative subsidence at time `t`, full grid.
"""
cumulative_subsidence(header::Header, t::Time) = cumulative_subsidence(header)(t)

"""
    cumulative_subsidence(header::Header, data::Data) -> Array{Length, D}

Per-cell cumulative subsidence at every write step of `data`, sliced to match
`data.slice`. The returned array has the same shape as `data.sediment_thickness`:
`(n_writes,)` for a column, `(n_pos, n_writes)` for a slice, `(n_x, n_y, n_writes)`
for a volume. Useful for plotters that need to broadcast subsidence against
other per-write quantities.
"""
function cumulative_subsidence(header::Header, data::Data{F, D}) where {F, D}
    cum = cumulative_subsidence(header)
    t_writes = header.axes.t[1:data.write_interval:end]
    out = similar(data.sediment_thickness, Length)
    for (k, t) in enumerate(t_writes)
        m = cum(t)                       # (nx, ny)
        selectdim(out, D, k) .= m[data.slice...]
    end
    return out
end

"""
    water_depth(header, data)

Compute the water depth function for the given data set.

The subsidence contribution is computed via `cumulative_subsidence(header, data)`,
which honours per-cell rate maps and modifiers when present.
"""
function water_depth(header::Header, data::Data{F, D}) where {F, D}
    na = [CartesianIndex()]
    sl = header.sea_level[1:data.write_interval:end]
    h0 = header.initial_topography[data.slice..., na]
    Δh = data.sediment_thickness
    Σ = cumulative_subsidence(header, data)
    return sl[repeated(na, D-1)...,:] .- h0 .- Δh .+ Σ
end

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
    set_attribute(out::T, name::String, value::Any)

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
    write_active_layer(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

Write the contents of the active layer to the output object.
"""
function write_active_layer end

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
function state_writer(input::Input, out) where {Input <: AbstractInput}
    output_sets = input.output
    grid_size = input.box.grid_size
    save_active_layer = hasfield(Input, :save_active_layer) ?
        input.save_active_layer : false

    return function (idx::Int, state::AbstractState)
        for (k, v) in output_sets
            if mod(idx - 1, v.write_interval) == 0
                write_sediment_thickness(
                    out, k, div(idx - 1, v.write_interval) + 1,
                    view(state.sediment_height, v.slice...))

                if save_active_layer
                    write_active_layer(
                        out, k, div(idx - 1, v.write_interval) + 1,
                        view(state.active_layer, :, v.slice...))
                end
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

    return function (idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, k, v) = ()
        function try_write(write::F, src, k, v) where {F}
            if (idx==1)
                write(out, k, 1,
                view(src, :, v.slice...))
            else
                write(out, k, div(idx - 2, v.write_interval) + 2,
                view(src, :, v.slice...))
            end
        end

        for (k, v) in input.output
            n_writes = div(input.time.steps, v.write_interval) + 1
            if div(idx-2, v.write_interval) + 2 <= n_writes
                try_write(write_production, frame.production, k, v)
                try_write(write_disintegration, frame.disintegration, k, v)
                try_write(write_deposition, frame.deposition, k, v)
            end
        end
    end
end

end
# ~/~ end
