# ~/~ begin <<docs/src/data-export.md#src/Export.jl>>[init]
module Export

export Data, DataSlice, DataColumn, Header, CSV, read_data, read_slice, read_column, data_export

using HDF5
import CSV: write as write_csv
using TOML

using Unitful
using DataFrames
using .Iterators: flatten

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

# ~/~ begin <<docs/src/data-export.md#export-specification>>[init]
abstract type ExportSpecification end

@kwdef struct CSV <: ExportSpecification
    output_files::IdDict{Symbol,String}
end

CSV(kwargs...) = CSV(IdDict(kwargs...))
# ~/~ end

@kwdef struct Axes
    x::Vector{Length}
    y::Vector{Length}
    t::Vector{Time}
end

const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}

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


function data_kind(gid::HDF5.Group)
	slice = parse_multi_slice(attrs(gid)["slice"])
	return data_kind(slice...)
end

function data_kind(fid::HDF5.File, group)
	group_name = string(group)
	if group_name == "input"
		return :metadata
	end
    return data_kind(fid[group_name])
end

function group_datasets(fid::HDF5.File)
	result = Dict{Symbol, Vector{String}}(
		:metadata => [],
		:volume => [],
		:slice => [],
		:column => [])

	for k in keys(fid)
		kind = data_kind(fid, k)
		push!(result[kind], k)
	end
	return result
end

function data_header(gid::HDF5.Group)
	slice = parse_multi_slice(attrs(gid)["slice"])
    kind = data_kind(slice...)
    write_interval = attrs(gid)["write_interval"]
    return DataHeader(
        slice=slice, kind=kind, write_interval=write_interval)
end

function read_header(fid)
    attrs = HDF5.attributes(fid["input"])

    axes = Axes(
        fid["input/x"][] * u"m",
        fid["input/y"][] * u"m",
        fid["input/t"][] * u"Myr")

    data_sets = Dict()
    for k in keys(fid)
        if k == "input"
            continue
        end
        data_sets[Symbol(k)] = data_header(fid[k])
    end

    return Header(
        attrs["tag"][],
        axes,
        attrs["delta_t"][] * u"Myr",
        attrs["time_steps"][],
        fid["input/initial_topography"][] * u"m",
        fid["input/sea_level"][] * u"m",
        attrs["subsidence_rate"][] * u"m/Myr",
        data_sets)
end

function read_data(::Type{Val{dim}}, gid::Union{HDF5.File, HDF5.Group}) where {dim}
	slice = parse_multi_slice(string(attrs(gid)["slice"]))
	write_interval = attrs(gid)["write_interval"]

	reduce(_) = (:)
	reduce(::Int) = 1

	Data{dim+1,dim}(
		slice, write_interval,
		gid["disintegration"][:, reduce.(slice)..., :] * u"m",
		gid["production"][:, reduce.(slice)..., :] * u"m",
		gid["deposition"][:, reduce.(slice)..., :] * u"m",
		gid["sediment_thickness"][reduce.(slice)..., :] * u"m")
end

function read_data(D::Type{Val{dim}}, filename::AbstractString, group) where {dim}
    h5open(filename) do fid
        header = read_header(fid)
		gid = fid[string(group)]
		data = read_data(D, gid)
        header, data
    end
end

read_volume(args...) = read_data(Val{3}, args...)
read_slice(args...) = read_data(Val{2}, args...)
read_column(args...) = read_data(Val{1}, args...)

time(header::Header, data::Data) = header.axes.t[1:data.write_interval:end]

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
    DataFrame((e.variable => df[!, e.variable] / unit(e.eltype)
               for e in eachrow(describe(df)))...)

"""
    write_unitful_csv(io::IO, df::DataFrame)

Write a CSV from a `DataFrame` with `Unitful` units. The units will be
represented in the CSV header, and stripped from the individual values.
"""
write_unitful_csv(io, df::DataFrame) =
    write_csv(io, ustrip(df), header=unitful_headers(df))
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[1]
Base.accumulate(f) = (args...; kwargs...) -> accumulate(f, args...; kwargs...)

"""
    age_depth_model(sediment_accumulation_curve::Vector)
    age_depth_model(sediment_accumulation_curve::DataFrame)

Compute the ADM from the SAC. Implemented as:

    reverse ∘ accumulate(min) ∘ reverse

The `DataFrame` version `select`s SAC columns, transformed into ADM.
"""
age_depth_model(sac::Vector{T}) where {T} = sac |> reverse |> accumulate(min) |> reverse
age_depth_model(sac_df::DataFrame) =
    let sac_cols = filter(contains("sac"), names(sac_df)),
        adm_cols = replace.(sac_cols, "sac" => "adm")

        select(sac_df, "timestep",
               (sac => age_depth_model => adm
                for (sac, adm) in zip(sac_cols, adm_cols))...)
    end
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[2]
"""
    stratigraphic_column(header::Header, column::DataColumn, facies::Int)

Compute the Stratigraphic Column for a given grid position `loc` and `facies` index.
Returns an `Array{Quantity, 2}` where the `Quantity` is in units of meters.
"""
function stratigraphic_column(header::Header, column::DataColumn, facies::Int)
	n_steps = size(column.production, 2)	
    delta = column.deposition[facies,:] .- column.disintegration[facies,:]

	for step in 1:n_steps
        debt = 0.0u"m"
        for pos in (step:-1:2)
            if delta[pos] > 0.0u"m"
                break
            end

            delta[pos-1] += delta[pos]
            delta[pos] = 0.0u"m"
        end
	end

	return delta
end
# ~/~ end
# ~/~ begin <<docs/src/data-export.md#export-function>>[3]
struct CSVExportTrait{S} end

"""
	data_export(spec::CSV, header::Header, data)

Exports `data` to CSV. Here, `data` should be a collection or iterable
of `DataColumn`.
"""
function data_export(spec::CSV, header::Header, data)
    for (key, filename) in spec.output_files
        if key == :metadata
            md = Dict(
                "global" => Dict(
                    "tag" => header.tag,
                    "subsidence_rate" => header.subsidence_rate,
                    "time_steps" => header.time_steps,
                    "delta_t" => header.Δt),
                "locations" => [Dict(
                    "label" => string(label),
                    "x" => header.axes.x[col.slice[1]],
                    "y" => header.axes.y[col.slice[2]],
                    "initial_topography" => header.initial_topography[col.slice...])
                    for (label, col) in pairs(data)],
                "files" => spec.output_files)
            open(filename, "w") do io
                TOML.print(io, md) do obj
                    if obj isa Quantity
                        [ustrip(obj), string(unit(obj))]
                    else
                        obj
                    end
                end
            end
            continue
        end
        open(filename, "w") do io
            time_df = DataFrame(
                :timestep => 0:header.time_steps,
                :time => header.axes.t)
            df = innerjoin(
                time_df,
                (data_export(CSVExportTrait{key}, header, column, label)
                 for (label, column) in pairs(data))...,
                on=:timestep)
            write_unitful_csv(io, df)
        end
    end
end

function data_export(::Type{CSVExportTrait{S}}, header::Header, data::DataColumn, label) where {S}
    error("Unknown CSV data export: `$(S)`")
end

function data_export(E::Type{CSVExportTrait{S}}, header::Header, columns) where {S}
    return innerjoin(
        (data_export(E, header, column, label)
         for (label, column) in pairs(columns))...,
        on=:timestep)
end

data_export(::Type{CSVExportTrait{:sediment_accumulation_curve}}, header::Header, data::DataColumn, label) =
    extract_sac(header, data, label)
data_export(::Type{CSVExportTrait{:stratigraphic_column}}, header::Header, data::DataColumn, label) =
    extract_sc(header, data, label)
data_export(::Type{CSVExportTrait{:water_depth}}, header::Header, data::DataColumn, label) =
    extract_wd(header, data, label)
data_export(::Type{CSVExportTrait{:age_depth_model}}, header::Header, data::DataColumn, label) =
    extract_sac(header, data, label) |> age_depth_model

"""
    extract_sac(header::Header, data::DataColumn)

Extract Sediment Accumumlation Curve (SAC) from the data. The SAC is directly
copied from `data.sediment_thickness`. Returns a `DataFrame` with `time` and
`sac_<n>` columns where `<n>` is in the range `1:length(grid_locations)`.
"""
function extract_sac(header::Header, data::DataColumn, label)
    DataFrame(
        "timestep" => 0:data.write_interval:header.time_steps, 
        "sac_$(label)" => data.sediment_thickness)
end

"""
    extract_sc(header::Header, data::DataColumn)

Extract Stratigraphic Column (SC) from the data. Returns a `DataFrame` with
`time` and `sc<n>` columns where `<n>` is in the range `1:length(grid_locations)`.
"""
function extract_sc(header::Header, data::DataColumn, label)
    n_facies = size(data.production)[1]
    DataFrame(
        "timestep" => data.write_interval:data.write_interval:header.time_steps, 
        ("sc_$(label)_f$(f)" => stratigraphic_column(header, data, f)
         for f in 1:n_facies)...)
end

"""
    extract_wd(header::Header, data::DataColumn)

Extract the water depth from the data. Returns a `DataFrame` with `time` and
`wd<n>` columns where `<n>` is in the range `1:length(grid_locations)`.
"""
function extract_wd(header::Header, data::DataColumn, label)
    na = [CartesianIndex()]
    t = header.axes.t[1:data.write_interval:end]
    sea_level = header.sea_level[1:data.write_interval:end]
    wd = header.subsidence_rate .* t .- 
        header.initial_topography[data.slice...] .- 
        data.sediment_thickness .+
        sea_level
    return DataFrame(
        "timestep" => 0:data.write_interval:header.time_steps, 
        "wd_$(label)" => wd)
end
# ~/~ end

end
# ~/~ end
