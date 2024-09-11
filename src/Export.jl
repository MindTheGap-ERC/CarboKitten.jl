# ~/~ begin <<docs/src/data-export.md#src/Export.jl>>[init]
module Export

using HDF5
using Unitful

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

struct Axes
	x::Vector{Length}
	y::Vector{Length}
	t::Vector{Time}
end

struct Header
	axes::Axes
	Δt::Time
	time_steps::Int
	bedrock_elevation::Matrix{Amount}
	sea_level::Vector{Length}
	subsidence_rate::Rate
end

struct Data
	disintegration::Array{Amount, 4}
	production::Array{Amount, 4}
	deposition::Array{Amount, 4}
	sediment_elevation::Array{Amount, 3}
end

struct DataSlice
	disintegration::Array{Amount, 3}
	production::Array{Amount, 3}
	deposition::Array{Amount, 3}
	sediment_elevation::Array{Amount, 2}
end

function read_header(fid)
	attrs = HDF5.attributes(fid["input"])

	axes = Axes(
		fid["input/x"][]*u"m",
		fid["input/y"][]*u"m",
		fid["input/t"][]*u"Myr")

	return Header(
		axes,
		attrs["delta_t"][]*u"Myr",
		attrs["time_steps"][],
		fid["input/bedrock_elevation"][]*u"m",
		fid["input/sea_level"][]*u"m",
		attrs["subsidence_rate"][]*u"m/Myr")
end

function read_data(filename)
	h5open(filename) do fid
		header = read_header(fid)
		data = Data(
			fid["disintegration"][]*u"m",
			fid["production"][]*u"m",
			fid["deposition"][]*u"m",
			fid["sediment_height"][]*u"m")
		header, data
	end
end

# ~/~ begin <<docs/src/data-export.md#extract-age-depth-model>>[init]
# ~/~ end

end
# ~/~ end