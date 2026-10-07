# ~/~ begin <<docs/src/output/data.md#src/Output/Storage.jl>>[init]
module Storage

using Unitful
using .Iterators: repeated

import ...CarboKitten: time_axis, box_axes, OutputSpec, AbstractOutput, AbstractInput, AbstractState
import ..Abstract: AbstractBundle, header

using ...Algorithms.StratigraphicColumn: stratigraphic_column!

export Header, Data, DataVolume, DataSlice, DataColumn, Pack, Axes, DataHeader
export data_kind, data_sets, data_volumes, data_slices, data_columns
export sediment_thickness, water_depth, surface_heights, stratigraphic_column

const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Slice2 = NTuple{2,Union{Int,Colon,UnitRange{Int}}}
const Amount = typeof(1.0u"m")
const Sediment = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")

# ~/~ begin <<docs/src/output/data.md#slice-helpers>>[init]
count_ints(::Int, args...) = 1 + count_ints(args...)
count_ints(_, args...) = count_ints(args...)
count_ints() = 0

reduce_slice(s::Tuple{Colon,Colon}, x, y) = (x, y)
reduce_slice(s::Tuple{Int,Colon}, y::Int) = (s[1], y)
reduce_slice(s::Tuple{Colon,Int}, x::Int) = (x, s[2])
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#slice-helpers>>[1]
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
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-header>>[init]
@kwdef struct Axes
    x::Vector{Length}
    y::Vector{Length}
    t::Vector{Time}
end
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-header>>[1]
@kwdef struct DataHeader
    kind::Symbol
    slice::Slice2
    write_interval::Int
end
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-header>>[2]
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
    data_sets::Dict{Symbol,DataHeader}
    attributes::Dict{String,Any} = Dict()
end

data_sets(header::Header) = header.data_sets
data_volumes(header::Header) = [s for (s, h) in header.data_sets if h.kind == :volume]
data_slices(header::Header) = [s for (s, h) in header.data_sets if h.kind == :slice]
data_columns(header::Header) = [s for (s, h) in header.data_sets if h.kind == :column]

data_sets(bundle::AbstractBundle) = bundle |> header |> data_sets
data_volumes(bundle::AbstractBundle) = bundle |> header |> data_volumes
data_slices(bundle::AbstractBundle) = bundle |> header |> data_slices
data_columns(bundle::AbstractBundle) = bundle |> header |> data_columns
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-data>>[init]
@kwdef struct Data{F,D}
    slice::Slice2
    write_interval::Int
    # Julia doesn't allow to say Array{Amount,D+1} here
    disintegration::Array{Amount,F}
    production::Array{Amount,F}
    deposition::Array{Amount,F}
    bathymetry::Array{Amount,D}
    active_layer::Union{Array{Amount,F}, Nothing} = nothing
    stratigraphic_column::Ref{Union{Array{Amount,F}, Nothing}} = nothing
end

const DataVolume = Data{4,3}
const DataSlice = Data{3,2}
const DataColumn = Data{2,1}

disintegration(v::Data) = v.disintegration
production(v::Data) = v.production
deposition(v::Data) = v.deposition
active_layer(v::Data) = v.active_layer
bathymetry(v::Data) = v.bathymetry
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-data>>[1]
Base.getindex(v::Data{F,D}, args...) where {F,D} =
    let k = count_ints(args...)
        Data{F - k,D - k}(
            reduce_slice(v.slice, args...),
            v.write_interval,
            v.disintegration[:, args..., :],
            v.production[:, args..., :],
            v.deposition[:, args..., :],
            v.bathymetry[args..., :],
            v.active_layer == nothing ? nothing : v.active_layer[:, args..., :],
            nothing)  # stratigraphic_column: reset so it is recomputed for the slice
    end
# ~/~ end
# ~/~ begin <<docs/src/output/data.md#data-pack>>[init]
struct Pack{F, D}
    header::Header
    data::Data{F, D}
end

bathymetry(v::Pack) = bathymetry(v.data)
production(v::Pack) = production(v.data)
deposition(v::Pack) = deposition(v.data)
active_layer(v::Pack) = active_layer(v.data)
disintegration(v::Pack) = disintegration(v.data)
# ~/~ end

"""
    stratigraphic_column(data)

Given a data set, compute the stratigraphic column. Result is memoised in
`data.stratigraphic_column` so repeated calls are free.
"""
function stratigraphic_column(data::Data{F, D}) where {F, D}
    if data.stratigraphic_column[] === nothing
        net_deposition = data.deposition .- data.disintegration
        for c in eachslice(net_deposition, dims=(1:D...,))
            stratigraphic_column!(c)
        end
        data.stratigraphic_column[] = net_deposition
    end
    return data.stratigraphic_column[]
end

stratigraphic_column(p::Pack{F, D}) where {F, D} = stratigraphic_column(p.data)

"""
    water_depth(header, data)

Compute the water depth function for the given data set.
"""
function water_depth(header::Header, data::Data{F, D}) where {F, D}
    sl = reshape(header.sea_level[1:data.write_interval:end], (repeated(1, D-1)..., :))
    return sl .- data.bathymetry
end

water_depth(p::Pack{F, D}) where {F, D} = water_depth(p.header, p.data)

"""
    sediment_thickness(data)

Compute the sediment thickness at each moment in the run by taking the cumulative
sum of the net deposition (deposition - disintegration) at each moment.
"""
function sediment_thickness(data::Data{F, D}) where {F, D}
    net_deposition = dropdims(sum(data.deposition .- data.disintegration, dims=1), dims=1)
    for c in eachslice(net_deposition, dims=(1:D-1...,))
        for i in 2:length(c)
            c[i] += c[i-1]
        end
    end
    return net_deposition
end

sediment_thickness(p::Pack{F, D}) where {F, D} = sediment_thickness(p.data)

"""
    surface_heights(header, data)

Compute the sediment surface height at every `(spatial..., time)` cell,
accounting for subsidence and net deposition. Returns an array of shape
`(spatial..., n_t+1)` where the first time entry is the initial topography
minus total subsidence and subsequent entries accumulate the preserved
sediment column.

Works with `DataColumn`, `DataSlice`, and `DataVolume`.
"""
function surface_heights(header::Header, data::Data{F, D}) where {F, D}
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    # Sum over the facies dimension (dim 1), yielding (spatial..., n_t)
    sc_sum = dropdims(sum(sc, dims=1), dims=1)
    # Cumulative sediment accumulation along the time axis (last dim)
    accumulated = cumsum(sc_sum, dims=ndims(sc_sum))
    n_t = size(sc_sum, ndims(sc_sum))
    h0 = initial_topography .- total_subsidence
    # Build result array: shape (spatial..., n_t+1)
    sz = (size(sc_sum)[1:end-1]..., n_t + 1)
    h = Array{eltype(h0)}(undef, sz...)
    selectdim(h, ndims(h), 1) .= h0
    selectdim(h, ndims(h), 2:n_t+1) .= h0 .+ accumulated
    return h
end

surface_heights(p::Pack{F, D}) where {F, D} = surface_heights(p.header, p.data)

end
# ~/~ end
