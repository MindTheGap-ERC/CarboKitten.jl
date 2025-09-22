# ~/~ begin <<docs/src/components/hdf5.md#src/Components/H5Writer.jl>>[init]
@compose module H5Writer
    using ..Common

    using HDF5

    import ...CarboKitten: run_model, Model
    import ...OutputData: AbstractOutputSpec, Frame

    @mixin Boxes, TimeIntegration, FaciesBase, WaterDepth

    # ~/~ begin <<docs/src/components/hdf5.md#hdf5-output-spec>>[init]
    @kwdef struct Input <: AbstractInput
        output = Dict(:full => OutputSpec((:,:), 1))
    end
    # ~/~ end

	axis_size(::Colon, a::Int) = a
	axis_size(::Int, _) = 1
	axis_size(r::AbstractRange{Int}, _) = length(r)

	slice_str(::Colon) = ":"
	slice_str(r::AbstractRange{Int}) = "$(r.start):$(r.stop)"
	slice_str(i::Int) = "$(i)"

    is_column(::Int, ::Int) = true
    is_column(_, _) = false

    function create_ck_group(fid, input::AbstractInput, name::Symbol, spec::AbstractOutputSpec)
        nf = n_facies(input)
        nw = div(input.time.steps, spec.write_interval)

        size = axis_size.(spec.slice, input.box.grid_size)

        grp = HDF5.create_group(fid, string(name))
        attrs = attributes(grp)
		attrs["slice"] = join(slice_str.(spec.slice),",")
        attrs["write_interval"] = spec.write_interval

        HDF5.create_dataset(grp, "production", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)
        HDF5.create_dataset(grp, "disintegration", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)
        HDF5.create_dataset(grp, "deposition", datatype(Float64),
            dataspace(nf, size..., nw),
            chunk=(nf, size..., 1), deflate=3)
        HDF5.create_dataset(grp, "sediment_thickness", datatype(Float64),
            dataspace(size..., nw + 1),
            chunk=(size..., 1), deflate=3)
    end

    function write_state(fid, input::AbstractInput, idx::Int, state::AbstractState)
        for (k, v) in input.output
			size = axis_size.(v.slice, input.box.grid_size)
            if mod(idx-1, v.write_interval) == 0
                fid[string(k)]["sediment_thickness"][:, :, div(idx-1, v.write_interval)+1] =
                    (is_column(v.slice...) ?
                        Float64[state.sediment_height[v.slice...] |> in_units_of(u"m");] :
                        reshape(state.sediment_height[v.slice...], size) |> in_units_of(u"m"))
            end
        end
    end

    function write_frame(fid, input::AbstractInput, idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, v) = ()
		n_f = n_facies(input)
        function try_write(tgt, src::AbstractArray, v)
			size = axis_size.(v.slice, input.box.grid_size)
            tgt[:, :, :, div(idx-1, v.write_interval) + 1] +=
				(reshape(src[:, v.slice...], (n_f, size...)) |> in_units_of(u"m"))
        end

        for (k, v) in input.output
            grp = fid[string(k)]
            try_write(grp["production"], frame.production, v)
            try_write(grp["disintegration"], frame.disintegration, v)
            try_write(grp["deposition"], frame.deposition, v)
        end
    end

    # ~/~ begin <<docs/src/components/hdf5.md#hdf5-run-model>>[init]
    """
        run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M

    Run a model and write output to HDF5. Here `M` should be a model, i.e. a
    module with `initial_state`, `step!` and `write_header` defined. Example:

        run_model(Model{ALCAP}, ALCAP.Example.INPUT, "example.h5")
    """
    function run_model(::Type{Model{M}}, input::AbstractInput, filename::AbstractString) where M
        state = M.initial_state(input)

        h5open(filename, "w") do fid
            create_group(fid, "input")
            M.write_header(fid, input)

            # create a group for every output item
            for (k, v) in input.output
                create_ck_group(fid, input, k, v)
            end
            write_state(fid, input, 1, state)

            run_model(Model{M}, input, state) do w, df
                # write_frame chooses to advance in a dataset
                # or just to increment on the current frame
                write_frame(fid, input, w, df)
                # write_state only writes one in every write_interval
                # and does no accumulation
                write_state(fid, input, w+1, state)
            end
        end

        return filename
    end
    # ~/~ end
end
# ~/~ end
