# ~/~ begin <<docs/src/components/hdf5.md#src/Components/H5Writer.jl>>[init]
@compose module H5Writer
    using ..Common
    using HDF5
    using ProgressLogging

    @mixin Boxes, TimeIntegration, FaciesBase, WaterDepth

    export run

    @kwdef struct DataFrame
        disintegration::Union{Array{Amount,3},Nothing} = nothing   # facies, x, y
        production::Union{Array{Amount,3},Nothing} = nothing
        deposition::Union{Array{Amount,3},Nothing} = nothing
    end

    Base.zeros(::Type{DataFrame}, input::AbstractInput) = DataFrame(
        disintegration=zeros(Amount, n_facies(input), input.box.grid_size...),
        production=zeros(Amount, n_facies(input), input.box.grid_size...),
        deposition=zeros(Amount, n_facies(input), input.box.grid_size...))

    function increment!(a::DataFrame, b::DataFrame)
        if !isnothing(b.disintegration)
            a.disintegration .+= b.disintegration
        end
        if !isnothing(b.production)
            a.production .+= b.production
        end
        if !isnothing(b.deposition)
            a.deposition .+= b.deposition
        end
    end

    function create_dataset(fid, input::AbstractInput)
        nf = n_facies(input)
        nw = n_writes(input)

        HDF5.create_dataset(fid, "production", datatype(Float64),
            dataspace(nf, input.box.grid_size..., nw),
            chunk=(nf, input.box.grid_size..., 1), deflate=3)
        HDF5.create_dataset(fid, "disintegration", datatype(Float64),
            dataspace(nf, input.box.grid_size..., nw),
            chunk=(nf, input.box.grid_size..., 1), deflate=3)
        HDF5.create_dataset(fid, "deposition", datatype(Float64),
            dataspace(nf, input.box.grid_size..., nw),
            chunk=(nf, input.box.grid_size..., 1), deflate=3)
        HDF5.create_dataset(fid, "sediment_height", datatype(Float64),
            dataspace(input.box.grid_size..., nw + 1),
            chunk=(input.box.grid_size..., 1), deflate=3)
    end

    function write_state(fid, idx::Int, state::AbstractState)
        fid["sediment_height"][:, :, idx] = state.sediment_height |> in_units_of(u"m")
    end

    function write_frame(fid, idx::Int, frame::DataFrame)
        fid["production"][:, :, :, idx] = frame.production |> in_units_of(u"m")
        fid["disintegration"][:, :, :, idx] = frame.disintegration |> in_units_of(u"m")
        fid["deposition"][:, :, :, idx] = frame.deposition |> in_units_of(u"m")
    end

    function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M
        state = M.initial_state(input)
        step! = M.step!(input)

        h5open(filename, "w") do fid
            create_group(fid, "input")
            M.write_header(fid, input)

            create_dataset(fid, input)
            write_state(fid, 1, state)
            @progress for w = 1:n_writes(input)
                df = zeros(DataFrame, input)
                for n = 1:input.time.write_interval
                    increment!(df, step!(state))
                end
                write_frame(fid, w, df)
                write_state(fid, w+1, state)
            end
        end

        return filename
    end
end
# ~/~ end
