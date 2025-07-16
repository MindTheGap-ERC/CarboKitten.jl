# ~/~ begin <<docs/src/memory-writer.md#src/OutputData.jl>>[init]
module OutputData

export Data, DataColumn, DataSlice, DataVolume, Slice2, Header, DataHeader, Axes
export parse_multi_slice, data_kind

using Unitful

const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}
const Amount = typeof(1.0u"m")
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
    initial_topography::Matrix{Amount}
    sea_level::Vector{Length}
    subsidence_rate::Rate
    data_sets::Dict{Symbol, DataHeader}
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

end
# ~/~ end
