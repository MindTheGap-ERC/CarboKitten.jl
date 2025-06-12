# Data export

We provide several ways to export reduced data from CarboKitten to CSV files that are easier to read, visualize and post-process for most people.

``` {.julia #export-specification}
abstract type ExportSpecification end

@kwdef struct CSV <: ExportSpecification
    output_files::IdDict{Symbol,String}
end

CSV(kwargs...) = CSV(IdDict(kwargs...))
```

```@example
using CarboKitten.Export: CSV

CSV(tuple.(10:20:70, 25),
  :sediment_accumulation_curve => "run_06_sac.csv",
  :age_depth_model             => "run_06_adm.csv",
  :stratigraphic_column        => "run_06_sc.csv",
  :water_depth                 => "run_06_wd.csv",
  :metadata                    => "run_06.toml")
```

There is a `data_export` function that can be overloaded with any `ExportSpcification`. When given a `CSV` specification, files are written as given.

- `:sediment_accumulation_curve`  (SAC) is another term for `sediment_height` elsewhere in the code.
- `:age_depth_model` (ADM) is a monotonic version of the SAC, relating depth to age.
- `:stratigraphic_column` amount of deposited material per facies per time step, corrected for disintegrated material. The cumulative sum of the SC should add up to the ADM.
- `:water_depth` the water depth as a function of time.
- `:metadata` some metadata, written as a TOML file.

## Tests

We have a test case with just three pixels.

1. Uniform production
2. Uniform production, top-hat disintegration, making the sediment accumulation non-monotonic
3. Linearly increasing production (not sure what this adds)

``` {.julia #export-test-case}
const AXES1 = Axes(
    x=[0.0, 1.0, 2.0] * u"m",
    y=[1.0] * u"m",
    t=(0.0:0.1:1.0) * u"Myr")

const HEADER1 = Header(
    tag="test",
    axes=AXES1,
    Δt=0.1u"Myr",
    time_steps=10,
    initial_topography=zeros(typeof(1.0u"m"), 3, 3),
    sea_level=zeros(typeof(1.0u"m"), 11),
    subsidence_rate=10u"m/Myr")

const PRODUCTION1 = reshape(
    hcat(ones(Amount, 10),
        ones(Amount, 10),
        cumsum(ones(Amount, 10)) / 5.5)',
    1, 3, 1, 10)

const DISINTEGRATION1 = reshape(
    hcat(zeros(Amount, 10),
        1:10 .|> (x -> x < 4 || x > 6 ? 0.0u"m" : 2.0u"m"),
        zeros(Amount, 10))',
    1, 3, 1, 10)

const ELEVATION1 = cat(
    [0.0, 0.0, 0.0]u"m",
    cumsum(PRODUCTION1 .- DISINTEGRATION1; dims=4)[1, :, :, :];
    dims=3)

const DATA1 = DataVolume(
    slice=(:,:),
    write_interval=1,
    disintegration=DISINTEGRATION1,
    production=PRODUCTION1,
    deposition=PRODUCTION1 .- DISINTEGRATION1,
    sediment_thickness=ELEVATION1)

const GRID_LOCATIONS1 = [(1, 1), (2, 1), (3, 1)]
```

### Sediment Accumulation

## Writing CSV files

The `CSV.jl` module lets us write a `DataFrame` to CSV, but doesn't work so well in combination with Units.

``` {.julia #export-function}
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
```

We may test that writing and reading the CSV back, gives the same result:

``` {.julia #export-test}
@testset "Hither and Dither" begin
    io = IOBuffer(UInt8[], read=true, write=true)
    data_export(CSVExportTrait{:sediment_accumulation_curve}, io, HEADER1, DATA1, GRID_LOCATIONS1)
    seek(io, 0)
    df = read_csv(io, DataFrame)
    rename!(df, (n => split(n)[1] for n in names(df))...)
    @test df.sac1 ≈ ELEVATION1[1, 1, :] / u"m"
    @test df.sac2 ≈ ELEVATION1[2, 1, :] / u"m"
    @test df.sac3 ≈ ELEVATION1[3, 1, :] / u"m"
end
```

## Age-depth model

``` {.julia #export-function}
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
```

We test that the constructed ADM is monotonic increasing in time:

``` {.julia #export-test}
@testset "ADM Monotonicity" begin
    sac = extract_sac(HEADER1, DATA1, GRID_LOCATIONS1)
    adm = sac |> age_depth_model

    @test sac.sac1 == adm.adm1
    @test sac.sac3 == adm.adm3
    @test sac.sac2 != adm.adm2

    @test all(adm.adm2[2:end] .- adm.adm2[1:end-1] .>= 0.0u"m")
end
```

## Stratigraphic Column

``` {.julia #export-function}
"""
    stratigraphic_column(header::Header, data::Data, loc::NTuple{2,Int}, facies::Int)

Compute the Stratigraphic Column for a given grid position `loc` and `facies` index.
Returns an `Array{Quantity, 2}` where the `Quantity` is in units of meters.
"""
function stratigraphic_column(header::Header, data::Data, loc::NTuple{2,Int}, facies::Int)
    dc = DataColumn(
        loc,
		data.write_interval,
        data.disintegration[:, loc..., :],
        data.production[:, loc..., :],
        data.deposition[:, loc..., :],
        data.sediment_thickness[loc..., :])
    return stratigraphic_column(header, dc, facies)
end

function stratigraphic_column(header::Header, data::DataColumn, facies::Int)
    n_times = length(header.axes.t) - 1
    sc = zeros(typeof(1.0u"m"), n_times)

    for ts = 1:n_times
        acc = data.deposition[facies, ts] - data.disintegration[facies, ts]
        if acc > 0.0u"m"
            sc[ts] = acc
            continue
        end
        ts_down = ts - 1
        while acc < 0.0u"m"
            ts_down < 1 && break
            if -acc < sc[ts_down]
                sc[ts_down] -= acc
                break
            else
                acc += sc[ts_down]
                sc[ts_down] = 0.0u"m"
            end
            ts_down -= 1
        end
    end

    sc
end
```

The stratigraphic column should sum to the age-depth model.

``` {.julia #export-test}
@testset "SC sum equals ADM" begin
    sac = extract_sac(HEADER1, DATA1, GRID_LOCATIONS1)
    adm = sac |> age_depth_model
    sc = extract_sc(HEADER1, DATA1, GRID_LOCATIONS1)
    @test [0.0u"m"; cumsum(sc.sc1_f1)] ≈ adm.adm1
    @test [0.0u"m"; cumsum(sc.sc2_f1)] ≈ adm.adm2
    @test [0.0u"m"; cumsum(sc.sc3_f1)] ≈ adm.adm3
end
```

## Dispatch

``` {.julia #export-function}
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
                    "number" => i,
                    "x" => header.axes.x[loc[1]],
                    "y" => header.axes.y[loc[2]],
                    "initial_topography" => header.initial_topography[loc...])
                                for (i, loc) in enumerate(spec.grid_locations)],
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

function data_export(::Type{CSVExportTrait{S}}, header::Header, data::DataColumn) where {S}
    error("Unknown CSV data export: `$(S)`")
end

data_export(::Type{CSVExportTrait{:sac}, header::Header, data::DataColumn, label) =
    extract_sac(header, data, label)
data_export(::Type{CSVExportTrait{:sc}, header::Header, data::DataColumn, label) =
    extract_sc(header, data, label)
data_export(::Type{CSVExportTrait{:wd}, header::Header, data::DataColumn, label) =
    extract_wd(header, data, label)
data_export(::Type{CSVExportTrait{:adm}, header::Header, data::DataColumn, label) =
    extract_sac(header, data, label) |> age_depth_model

"""
    extract_sac(header::Header, data::DataColumn)

Extract Sediment Accumumlation Curve (SAC) from the data. The SAC is directly
copied from `data.sediment_thickness`. Returns a `DataFrame` with `time` and
`sac_<n>` columns where `<n>` is in the range `1:length(grid_locations)`.
"""
function extract_sac(header::Header, data::DataColumn, label)
    DataFrame(
        "timestep" => 0:header.time_steps,
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
        "timestep" => 1:header.time_steps, 
        ("sc$(i)_f$(f)" => stratigraphic_column(header, data, f)
         for f in 1:n_facies)...)
end

"""
    extract_wd(header::Header, data::DataColumn)

Extract the water depth from the data. Returns a `DataFrame` with `time` and
`wd<n>` columns where `<n>` is in the range `1:length(grid_locations)`.
"""
function extract_wd(header::Header, data::DataColumn, label)
    na = [CartesianIndex()]
    t = header.time[1:data.write_interval:end]
    sea_level = header.sea_level[1:data.write_interval:end]
    wd = header.subsidence_rate .* t .- 
        header.initial_topography[data.slice...] .- 
        data.sediment_thickness[data.slice...] .+
        sea_level
    return DataFrame(
        "timestep" => 0:header.time_steps, 
        "wd_$(label)" => wd)
end
```

``` {.julia file=src/Export.jl}
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

<<export-specification>>

@kwdef struct Axes
    x::Vector{Length}
    y::Vector{Length}
    t::Vector{Time}
end

@kwdef struct Header
    tag::String
    axes::Axes
    Δt::Time
    time_steps::Int
    initial_topography::Matrix{Amount}
    sea_level::Vector{Length}
    subsidence_rate::Rate
end

const Slice2 = NTuple{2, Union{Int, Colon, UnitRange{Int}}}

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

function data_kind(fid::HDF5.File, group)
	group_name = string(group)
	if group_name == "input"
		return :metadata
	end
	gid = fid[group_name]
	slice = parse_multi_slice(attrs(gid)["slice"])
	return data_kind(slice...)
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

function read_header(fid)
    attrs = HDF5.attributes(fid["input"])

    axes = Axes(
        fid["input/x"][] * u"m",
        fid["input/y"][] * u"m",
        fid["input/t"][] * u"Myr")

    return Header(
        attrs["tag"][],
        axes,
        attrs["delta_t"][] * u"Myr",
        attrs["time_steps"][],
        fid["input/initial_topography"][] * u"m",
        fid["input/sea_level"][] * u"m",
        attrs["subsidence_rate"][] * u"m/Myr")
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

<<export-function>>

end
```

``` {.julia file=test/ExportSpec.jl}
using CarboKitten
using CarboKitten.Export: Axes, Header, Data, data_export, CSVExportTrait,
    age_depth_model, extract_sac, extract_sc, CSV, read_data, extract_sac, extract_wd
using CSV: read as read_csv
using TOML
using DataFrames
using Unitful

const Amount = typeof(1.0u"m")

<<export-test-case>>

@testset "Data Export" begin
    <<export-test>>

    @testset "Write to folder" begin
        mktempdir() do path
            spec = CSV(GRID_LOCATIONS1,
                :sediment_accumulation_curve => joinpath(path, "sac.csv"),
                :age_depth_model => joinpath(path, "adm.csv"),
                :stratigraphic_column => joinpath(path, "sc.csv"),
                :water_depth => joinpath(path, "wd.csv"),
                :metadata => joinpath(path, "metadata.toml"))
            data_export(spec, HEADER1, DATA1)
            for f in values(spec.output_files)
                @test isfile(f)
            end

            metadata = TOML.parsefile(spec.output_files[:metadata])
            @test IdDict(Symbol(k) => v for (k, v) in metadata["files"]) == spec.output_files
            @test length(metadata["locations"]) == 3
            adm = read_csv(spec.output_files[:age_depth_model], DataFrame)
            rename!(adm, (n => split(n)[1] for n in names(adm))...)
            @test adm == ustrip(extract_sac(HEADER1, DATA1, GRID_LOCATIONS1) |> age_depth_model)
        end
    end

    @testset "Waterdepth signs" begin
        BS92_TEST_INPUT = BS92.Input(
            tag = "single pixel model",
            box = Box{Periodic{2}}(grid_size=(1, 1), phys_scale=600.0u"m"),
            time = TimeProperties(
              Δt = 10.0u"yr",
              steps = 8000,
              write_interval = 100),
            sea_level = t -> 10.0u"m" * sin(2π * t / 20u"kyr"),
            initial_topography = (_, _) -> - 50.0u"m",
            subsidence_rate = 0.001u"m/yr",
            insolation = 400.0u"W/m^2",
            facies = [BS92.Facies(
              maximum_growth_rate = 0.005u"m/yr",
              saturation_intensity = 50.0u"W/m^2",
              extinction_coefficient = 0.05u"m^-1"
            )])

        mktempdir() do path
            run_model(Model{BS92}, BS92_TEST_INPUT, joinpath(path, "run.h5"))
            header, data = read_data(joinpath(path, "run.h5"))
            wd = extract_wd(header, data, [(1, 1)])
            sac = extract_sac(header, data, [(1, 1)])
            submerged = wd.wd1 .> -1.0u"m"
            growing = (sac.sac1[2:end] .- sac.sac1[1:end-1]) .> 0.5u"m"
            @test all(growing .&& (submerged[1:end-1] .|| submerged[2:end]) .|| .!growing)
        end
    end
end
```
