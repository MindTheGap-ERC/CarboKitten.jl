# ~/~ begin <<docs/src/output/h5reader.md#src/Output/H5Reader.jl>>[init]
module H5Reader

using HDF5
import ..Abstract: AbstractBundle
import ..Storage: data_kind, Header, DataHeader, Data, Pack, data_sets, data_volumes, data_slices, data_columns

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

    grid_size = (length(axes.x), length(axes.y))
    n_facies = attrs["n_facies"][]

    return Header(
        tag = attrs["tag"][],
        axes = axes,
        Δt = attrs["delta_t"][] * u"Myr",
        time_steps = attrs["time_steps"][],
        grid_size = grid_size,
        n_facies = n_facies,
        initial_topography = fid["input/initial_topography"][] * u"m",
        sea_level = fid["input/sea_level"][] * u"m",
        subsidence_rate = attrs["subsidence_rate"][] * u"m/Myr",
        data_sets = data_sets)
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
		gid["bathymetry"][reduce.(slice)..., :] * u"m",
		"active_layer" in keys(gid) ?
		    gid["active_layer"][:, reduce.(slice)..., :] * u"m" :
			nothing, nothing)
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

struct H5Bundle <: AbstractBundle
    fid::HDF5.File
    header::Header
end

Base.close(bundle::H5Bundle) = close(bundle.fid)
header(bundle::H5Bundle) = bundle.header

function load(filename::AbstractString)
    fid = h5open(filename)
    header = read_header(fid)
    bundle = H5Bundle(fid, header)
    return bundle
end

function load_group(::Val{D}, bundle::H5Bundle, group) where {D}
    gid = bundle.fid[string(group)]
	data = read_data(D, gid)
    return Pack(bundle.header, data)
end

load_volume(bundle::H5Bundle, group::Symbol) = load_group(Val{3}, bundle, group)
load_volume(filename::AbstractString, group::Symbol) = load_volume(load(filename), group)
load_slice(bundle::H5Bundle, group::Symbol) = load_group(Val{2}, bundle, group)
load_slice(filename::AbstractString, group::Symbol) = load_slice(load(filename), group)
load_column(bundle::H5Bundle, group::Symbol) = load_group(Val{1}, bundle, group)
load_column(filename::AbstractString, group::Symbol) = load_column(load(filename), group)

end
# ~/~ end
