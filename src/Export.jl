# ~/~ begin <<docs/src/data-export.md#src/Export.jl>>[init]
module Export

using HDF5
using Unitful
using DataFrames
import CSV: write as write_csv

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

# ~/~ begin <<docs/src/data-export.md#export-specification>>[init]
abstract type ExportSpecification end

@kwdef struct CSV <: ExportSpecification
    grid_locations::Vector{NTuple{2,Int}}
    output_files::IdDict{Symbol,String}
end

CSV(; grid_locations, kwargs...) = CSV(grid_locations, kwargs)
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
    disintegration::Array{Amount,4}
    production::Array{Amount,4}
    deposition::Array{Amount,4}
    sediment_elevation::Array{Amount,3}
end

struct DataSlice
    disintegration::Array{Amount,3}
    production::Array{Amount,3}
    deposition::Array{Amount,3}
    sediment_elevation::Array{Amount,2}
end

function read_header(fid)
    attrs = HDF5.attributes(fid["input"])

    axes = Axes(
        fid["input/x"][] * u"m",
        fid["input/y"][] * u"m",
        fid["input/t"][] * u"Myr")

    return Header(
        axes,
        attrs["delta_t"][] * u"Myr",
        attrs["time_steps"][],
        fid["input/bedrock_elevation"][] * u"m",
        fid["input/sea_level"][] * u"m",
        attrs["subsidence_rate"][] * u"m/Myr")
end

function read_data(filename)
    h5open(filename) do fid
        header = read_header(fid)
        data = Data(
            fid["disintegration"][] * u"m",
            fid["production"][] * u"m",
            fid["deposition"][] * u"m",
            fid["sediment_height"][] * u"m")
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
    DataFrame((e.variable => df[!, e.variable] / unit(e.eltype)
               for e in eachrow(describe(df)))...)

"""
    write_unitful_csv(io::IO, df::DataFrame)

Write a CSV from a `DataFrame` with `Unitful` units. The units will be
represented in the CSV header, and stripped from the individual values.
"""
write_unitful_csv(io, df::DataFrame) =
    write_csv(io, ustrip(df), header=unitful_headers(df))

Base.accumulate(f) = (args...; kwargs...) -> accumulate(f, args...; kwargs...)

"""
    age_depth_curve(sediment_accumulation_curve)

Compute the ADM from the SAC. Implemented as:

    reverse ∘ accumulate(min) ∘ reverse
"""
age_depth_model(sac::Vector{T}) where {T} = sac |> reverse |> accumulate(min) |> reverse
age_depth_model(sac_df::DataFrame) =
    let sac_cols = filter(contains("sac"), names(sac_df)),
        adm_cols = replace.(sac_cols, "sac" => "adm")

        select(sac_df, "time", (sac => age_depth_model => adm
                                for (sac, adm) in zip(sac_cols, adm_cols))...)
    end

struct CSVExportTrait{S} end

function data_export(spec::T, filepath::String) where {T<:ExportSpecification}
    data_export(spec, read_data(filepath)...)
end

function data_export(spec::CSV, header::Header, data::Data)
    for (key, filename) in spec.output_files
        open(filename, "w") do io
            data_export(CSVExportTrait{key}, io, header, data, spec.grid_locations)
        end
    end
end

function data_export(::Type{CSVExportTrait{S}}, args...) where {S}
    error("Unknown CSV data export: `$(S)`")
end

function data_export(::Type{CSVExportTrait{:sac}}, io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    sac = extract_sac(header, data, grid_locations)
    write_unitful_csv(io, sac)
end

function data_export(::Type{CSVExportTrait{:adm}}, io::IO, header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    adm = extract_sac(header, data, grid_locations) |> age_depth_model
    write_unitful_csv(io, adm)
end

function extract_sac(header::Header, data::Data, grid_locations::Vector{NTuple{2,Int}})
    n_cols = length(grid_locations)
    DataFrame(:time => header.axes.t[2:end],
        (Symbol("sac$(i)") => data.sediment_elevation[grid_locations[i]..., :]
         for i = 1:n_cols)...)
end
# ~/~ end

end
# ~/~ end
