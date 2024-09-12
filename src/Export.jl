# ~/~ begin <<docs/src/data-export.md#src/Export.jl>>[init]
module Export

using HDF5
using Unitful
using DataFrames

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

# ~/~ begin <<docs/src/data-export.md#export-specification>>[init]
abstract type ExportSpecification end

struct CSV <: ExportSpecification
    grid_locations::Vector{NTuple{2,Int}}
    output_files::IdDict{Symbol,String}
end
# ~/~ end

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

# ~/~ begin <<docs/src/data-export.md#export-function>>[init]
"""
    unitful_headers(df::DataFrame)

Gets a string representation for all column names including their unit.
Returns a `Vector{String}`.
"""
unitful_headers(df::DataFrame) =
	["$(e.variable) [$(unit(e.eltype))]" for e in eachrow(describe(df))]

"""
    ustrip(df::DataFrame)

Strip units from a `DataFrame`. Returns a new `DataFrame`.
"""
Unitful.ustrip(df::DataFrame) =
	DataFrame((e.variable => df[!,e.variable] / unit(e.eltype)
		for e in eachrow(describe(df)))...)

"""
    write_unitful_csv(io::IO, df::DataFrame)

Write a CSV from a `DataFrame` with `Unitful` units. The units will be
represented in the CSV header, and stripped from the individual values.
"""
write_unitful_csv(io::IO, df::DataFrame) =
	CSV.write(io, ustrip(df), header=unitful_headers(df))

struct CSVExportTrait{S} end

function data_export(spec::CSV, filepath::String)
    header, data = read_data(file_path)
    for (key, filename) in spec.output_files
        open(filename, "w") do io
            data_export(CSVEXportTrait{key}, io, header, data, spec.grid_locations)
        end
    end
end

function data_export(::Type{CSVExportTrait{S}}, args...) where {S}
    error("Unknown CSV data export: `$(S)`")
end

function data_export(::Type{CSVExportTrait{:sac}}, io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2, Int}})
    df = extract_sac(header, data, grid_locations)
    write_unitful_csv(io, df)
end

function extract_sac(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    n_cols = length(grid_locations)
    # column_names = ["time (Myr)", ("sac$(i) (m)" for i in 1:n_cols)]
    DataFrame(:time => header.axes.t[2:end],
        (Symbol("sac$(i)") => data.sediment_elevation[grid_locations[i]...,:]
            for i = 1:n_cols)...)
end
# ~/~ end

end
# ~/~ end
