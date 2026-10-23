IO Interface
============

Writers
-------

The following defines the abstract writer interface. There are currently two implementations: `H5Writer` and `MemoryOutput`.

``` {.julia #abstract-writer}
const Sediment = typeof(1.0u"m")

@kwdef struct Frame
    disintegration::Union{Array{Sediment,3},Nothing} = nothing   # facies, x, y
    production::Union{Array{Sediment,3},Nothing} = nothing
    deposition::Union{Array{Sediment,3},Nothing} = nothing
end

"""
    new_output(::Type{T}, input)

Create a new output object of type `T`, given `input`.
"""
function new_output end

"""
    add_data_set(out::T, name::Symbol, spec::OutputSpec)

Add a data set to the output object.
"""
function add_data_set end

"""
    set_attribute(out::T, name::String, value::Any)

Set an attribute in the output object.
"""
function set_attribute end

"""
    write_bathymetry(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

Write the bathymetry to the output object. The `idx` should be corrected for
write interval. That is, `idx` should range from `1` to `n_writes` for the named
data set. This function should be implemented for 0, 1, and 2 dimensional
arrays, corresponding to writing column, slice or volume data.

If your output object type doesn't conform to the standard CarboKitten data
layout, you may choose to not implement this function and implement
`state_writer` and `frame_writer` instead. The same goes for `write_production`,
`write_disintegration` and `write_deposition`.
"""
function write_bathymetry end

"""
    write_active_layer(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

Write the contents of the active layer to the output object.
"""
function write_active_layer end

"""
    write_production(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_production end

"""
    write_disintegration(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_disintegration end

"""
    write_deposition(out::T, name::Symbol, idx::Int, data::AbstractArray{Amount, dim}) where {T, dim}

See `write_sediment_thickness`. Should accept 1, 2, and 3 dimensional arrays, corresponding to
writing column, slice or volume data. (first axis is facies, then x and y)
"""
function write_deposition end

"""
    state_writer(input::AbstractInput, out::T)

Returns a `function (idx::Int, state::AbstractState)`.

Write the state of the simulation to the output object. It is the responsibility
of the implementation to choose to write or not based on the set write interval.

The default implementation writes the state for all output data sets, and calls
`write_sediment_thickness`.
"""
function state_writer(input::Input, out) where {Input <: AbstractInput}
    output_sets = input.output
    grid_size = input.box.grid_size
    save_active_layer = hasfield(Input, :save_active_layer) ?
        input.save_active_layer : false

    return function (idx::Int, state::AbstractState)
        for (k, v) in output_sets
            if mod(idx - 1, v.write_interval) == 0
                write_bathymetry(
                    out, k, div(idx - 1, v.write_interval) + 1,
                    view(state.bathymetry, v.slice...))

                if save_active_layer
                    write_active_layer(
                        out, k, div(idx - 1, v.write_interval) + 1,
                        view(state.active_layer, :, v.slice...))
                end
            end
        end
    end
end

"""
    frame_writer(input::AbstractInput, out::T)

Returns a `function (idx::Int, state::AbstractState)`.

Write the state of the simulation to the output object. It is the responsibility
of the implementation to choose to write or not based on the set write interval.

The default implementation writes the state for all output data sets, and calls
`write_sediment_thickness`.
"""
function frame_writer(input::AbstractInput, out)
    n_f = length(input.facies)
    grid_size = input.box.grid_size

    return function (idx::Int, frame::Frame)
        try_write(tgt, ::Nothing, k, v) = ()
        function try_write(write::F, src, k, v) where {F}
            if (idx==1)
                write(out, k, 1,
                view(src, :, v.slice...))
            else
                write(out, k, div(idx - 2, v.write_interval) + 2,
                view(src, :, v.slice...))
            end
        end

        for (k, v) in input.output
            n_writes = div(input.time.steps, v.write_interval) + 1
            if div(idx-2, v.write_interval) + 2 <= n_writes
                try_write(write_production, frame.production, k, v)
                try_write(write_disintegration, frame.disintegration, k, v)
                try_write(write_deposition, frame.deposition, k, v)
            end
        end
    end
end
```

Readers
-------

The following defines the abstract reader interface. An implementation of `AbstractBundle` should also implement `Base.close`.

``` {.julia #abstract-reader}
"""
    load(filename::AbstractString)

Load data from `filename`. Returns a `H5Bundle`.

    load(output::AbstractOutput)

Load data from `output`. Returns an `AbstractBundle`.

    load(f::Function, args...)

Load data using `f` and `args`, and close the bundle after use. Use this
with a `do` block.
"""
function load end

function load(f::Function, args...)
    bundle = load(args...)
    try
        result = f(bundle)
        return result
    finally
        close(bundle)
    end
end

"""
    load_volume(bundle, sym)

Load a volume from `bundle`.
"""
function load_volume end

"""
    load_slice(bundle, sym)

Load a slice from `bundle`.
"""
function load_slice end

"""
    load_column(bundle, sym)

Load a column from `bundle`.
"""
function load_column end

abstract type AbstractBundle end

"""
    header(bundle)

Get the header of `bundle`.
"""
function header end
```

Module
------

``` {.julia file=src/Output/Abstract.jl}
module Abstract

import ...CarboKitten: set_attribute  # TODO: get rid of this

export Frame, new_output, add_data_set, set_attribute, state_writer, frame_writer
export write_bathymetry, write_active_layer, write_production, write_deposition, write_disintegration
export AbstractBundle, load, load_volume, load_slice, load_column, header

using Unitful
using ...CarboKitten: AbstractInput, AbstractState

<<abstract-writer>>
<<abstract-reader>>

end
```

Run model
---------

On top of this we have defined a `run_model` method that writes output in some form.

``` {.julia #run-model-output}
"""
    run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where M

Run a model and save the output to `output`.
"""
function run_model(::Type{Model{M}}, input::AbstractInput, output::AbstractOutput) where {M}
    M.write_header(input, output)

    state = M.initial_state(input)
    write_state = state_writer(input, output)
    write_frame = frame_writer(input, output)

    # create a group for every output item
    for (k, v) in input.output
        add_data_set(output, k, v)
    end
    write_state(1, state)
    # also write any initial sediment to output
    write_frame(1, M.initial_frame(input))

    run_model(Model{M}, input, state) do w, df
        # write_frame chooses to advance in a dataset
        # or just to increment on the current frame
        write_frame(w + 1, df)
        # write_state only writes one in every write_interval
        # and does no accumulation
        write_state(w + 1, state)
    end

    return output
end
```

``` {.julia file=src/Output/RunModel.jl}
module RunModel

import ...CarboKitten: run_model, Model
using ...CarboKitten: AbstractInput, AbstractOutput
using ..Abstract

<<run-model-output>>

end
```

## Tests

We had a recurring bug where unconformities are plotted wrong in cases of extreme erosion. It turned out, the function that computes water depth contained a subtle bug. The following test runs a model with extreme erosion and verifies that the bathymetry and water depth is computed correctly.

``` {.julia file=test/Output/WaterDepthSpec.jl}
using CarboKitten
using CarboKitten.Models: WithoutCA as M

function test_model_input()
    facies = [
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:euphotic] * ProductionBoost(factor=0.25),
            transport_coefficient=50.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:oligophotic] * ProductionBoost(factor=0.25),
            transport_coefficient=25.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:aphotic] * ProductionBoost(factor=0.25),
            transport_coefficient=12.5u"m/Myr"),
    ]

    M.Input(
        box = CarboKitten.Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
        time = TimeProperties(Δt=0.0002u"Myr", steps=5000),
        facies = facies,
        initial_topography = (x, y) -> -x / 300.0,

        sea_level = t -> 5.5u"m" * sin(2π * t / 0.18u"Myr"),
        subsidence_rate = 38.0u"m/Myr",
        disintegration_rate = 80.0u"m/Myr",
        lithification_time = 100.0u"yr",
        insolation = 400.0u"W/m^2",
    )
end

function test_model()
    input = test_model_input()
    output = MemoryOutput(input)
    run_model(Model{M}, input, output)

end

@testset "Output.Abstract.water_depth" begin
    na = [CartesianIndex()]
    output = test_model()
    section = output.data_volumes[:full][:, 1]
    wd = water_depth(output.header, section)
    sc = stratigraphic_column(section)
    st_preserved = dropdims(cumsum(sum(sc; dims=1); dims=3), dims=1)
    st = dropdims(
        cumsum(sum(section.deposition .- section.disintegration; dims=1), dims=3),
        dims=1)

    @test st_preserved[:, end] ≈ st[:, end]

    initial_bathymetry = output.header.initial_topography[:, na]
    subsidence = output.header.subsidence_rate .* output.header.axes.t[na, :]
    reconstructed_bathymetry = initial_bathymetry .+ st .- subsidence

    @test reconstructed_bathymetry ≈ section.bathymetry

    sea_level = output.header.sea_level[na, :]
    reconstructed_wd = sea_level .- reconstructed_bathymetry

    @test reconstructed_wd ≈ wd
end
```
