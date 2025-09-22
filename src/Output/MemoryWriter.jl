# ~/~ begin <<docs/src/memory-writer.md#src/Output/MemoryWriter.jl>>[init]
module MemoryWriter

using ..OutputData
import ..OutputData:
    new_output, add_data_set, set_attribute, write_sediment_thickness,
    write_production, write_disintegration, write_deposition, AbstractOutputSpec

using ..Components.Common
using ..Components.WaterDepth: initial_topography
using ..CarboKitten: time_axis, box_axes

struct MemoryOutput <: AbstractOutput
    header::Header
    data_volumes::Dict{Symbol, DataVolume}
    data_slices::Dict{Symbol, DataSlice}
    data_columns::Dict{Symbol, DataColumn}
end

MemoryOutput(input::AbstractInput) = new_output(MemoryOutput, input)

function new_output(::Type{MemoryOutput}, input::AbstractInput)
    t_axis = time_axis(input.time)
    x_axis, y_axis = box_axes(input.box)
    axes = Axes(x = x_axis, y = y_axis, t = t_axis)
    h0 = initial_topography(input)
    sl = input.sea_level.(t_axis)

    header = Header(
        tag = input.tag,
        axes = axes,
        Δt = input.time.Δt,
        time_steps = input.time.steps,
        grid_size = input.box.grid_size,
        n_facies = length(input.facies),
        initial_topography = h0,
        sea_level = sl,
        subsidence_rate = input.subsidence_rate,
        data_sets = Dict(),
        attributes = Dict())

    return MemoryOutput(header, Dict(), Dict(), Dict())
end

axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)

function add_data_set(out::MemoryOutput, label::Symbol, spec::AbstractOutputSpec)
    h = DataHeader(data_kind(spec), spec.slice, spec.write_interval)
    out.header.data_sets[label] = h

    slice = spec.slice
    write_interval = spec.write_interval

    full_size = out.header.grid_size
    n_steps = div(out.header.time_steps, write_interval)
    n_facies = out.header.n_facies

    if h.kind == :volume
        size = axis_size.(slice, full_size)
        out.data_volumes[label] = DataVolume(
            slice, write_interval,
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, n_facies, size..., n_steps),
            zeros(Amount, size..., n_steps + 1))
    elseif h.kind == :slice
        size = axis_size.(slice, full_size)
        slice_size = size[1] == 1 ? size[2] : size[1]
        out.data_slices[label] = DataSlice(
            slice, write_interval,
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, n_facies, slice_size, n_steps),
            zeros(Amount, slice_size, n_steps + 1))
    elseif h.kind == :column
        out.data_columns[label] = DataColumn(
            slice, write_interval,
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_facies, n_steps),
            zeros(Amount, n_steps + 1))
    end
end

function set_attribute(out::MemoryOutput, name::Symbol, value::Any)
    out.header.attributes[name] = value
end

write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 0}) =
    out.data_columns[label].sediment_thickness[idx] = data[]
write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_slices[label].sediment_thickness[:, idx] .= data
write_sediment_thickness(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_volumes[label].sediment_thickness[:, :, idx] .= data

write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].production[:, idx] .+= data
write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].production[:, :, idx] .+= data
write_production(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].production[:, :, :, idx] .+= data

write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].disintegration[:, idx] .+= data
write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].disintegration[:, :, idx] .+= data
write_disintegration(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].disintegration[:, :, :, idx] .+= data

write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 1}) =
    out.data_columns[label].deposition[:, idx] .+= data
write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 2}) =
    out.data_slices[label].deposition[:, :, idx] .+= data
write_deposition(out::MemoryOutput, label::Symbol, idx::Int, data::AbstractArray{Amount, 3}) =
    out.data_volumes[label].deposition[:, :, :, idx] .+= data

end
# ~/~ end