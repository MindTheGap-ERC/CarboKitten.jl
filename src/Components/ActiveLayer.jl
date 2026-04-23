# ~/~ begin <<docs/src/active-layer-transport.md#src/Components/ActiveLayer.jl>>[init]
@compose module ActiveLayer
@mixin WaterDepth, FaciesBase, SedimentBuffer

export disintegrator, transporter, cementation_factor

using ..Common
using ..Production
using CarboKitten.Transport.Advection: transport, advection_coef!, transport_dC!, max_dt
using CarboKitten.Transport.Solvers: runge_kutta_4, forward_euler
using Unitful
using GeometryBasics



@kwdef struct Facies <: AbstractFacies
    diffusion_coefficient::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    wave_velocity::Function = _ -> (
        Vec2(0.0u"m/Myr", 0.0u"m/Myr"),
        Vec2(0.0u"1/Myr", 0.0u"1/Myr")
    )

    depositional_porosity::Float64 = 0.55
    compaction_curve::Function = _ -> 0.0
end

@kwdef mutable struct State <: AbstractState
            active_layer::Array{Amount,3}
    compaction::Production.CompactionState = Production.CompactionState(0, 0, 1)

    layer_facies_hist::Vector{Matrix{Int16}} = Matrix{Int16}[]
    layer_thickness_hist::Vector{Matrix{Float64}} = Matrix{Float64}[]
    layer_wdepth_hist::Vector{Matrix{Float32}} = Matrix{Float32}[]
    layer_energy_hist::Vector{Matrix{Float32}} = Matrix{Float32}[]
end

@kwdef struct Input <: AbstractInput
    intertidal_zone::Height = 0.0u"m"
    disintegration_rate::Rate = 50.0u"m/Myr"
    disintegration_transfer::Function = x -> x
    cementation_time::Union{typeof(1.0u"Myr"),Nothing} = nothing
    transport_solver = Val{:RK4}
    transport_substeps = :adaptive
    wave_velocity::Union{Nothing, Function} = nothing
end

courant_max(::Type{Val{:RK4}}) = 2.0
courant_max(::Type{Val{:forward_euler}}) = 1.0

transport_solver(f, _) = f
transport_solver(::Type{Val{:RK4}}, box) = runge_kutta_4(typeof(1.0u"m"), box)
transport_solver(::Type{Val{:forward_euler}}, _) = forward_euler

function cementation_factor(input::AbstractInput)
    if input.cementation_time === nothing
        return 1.0
    else
        return 1.0 - exp(input.time.Δt * log(1 / 2) / input.cementation_time)
    end
end


# ------------------------------------------------------------------
# Adaptive transporter (now supports spatial wave fields)
# ------------------------------------------------------------------

function adaptive_transporter(input)

    solver = transport_solver(input.transport_solver, input.box)

    w = water_depth(input)
    box = input.box
    Δt = input.time.Δt
    fs = input.facies

    adv = Matrix{Vec2{Rate}}(undef, box.grid_size...)
    rct = Matrix{typeof(1.0u"1/Myr")}(undef, box.grid_size...)
    dC = Matrix{Rate}(undef, box.grid_size...)

    cm = courant_max(input.transport_solver)
    iz = input.intertidal_zone

    return function (state)

        wd = w(state)
        wd .+= iz

        C = state.active_layer

        for (i, f) in pairs(fs)

            # --------------------------------------------------
            # Spatial wave fields if present (ALCAP only)
            # --------------------------------------------------
            if hasproperty(state, :wave_velocity_field) &&
               hasproperty(state, :wave_shear_field) &&
               input.wave_velocity === nothing

                advection_coef!(box,
                                f.diffusion_coefficient,
                                (state.wave_velocity_field,
                                 state.wave_shear_field),
                                wd, adv, rct)

            else
                wv = isnothing(input.wave_velocity) ?
                     f.wave_velocity :
                     input.wave_velocity

                advection_coef!(box,
                                f.diffusion_coefficient,
                                wv, wd, adv, rct)
            end

            m = max_dt(adv, box.phys_scale, cm)
            steps = ceil(Int, Δt / m)
            subdt = Δt / steps

            for _ in 1:steps
                solver(
                    (C, _) -> transport_dC!(box, adv, rct, C, dC),
                    view(C, i, :, :),
                    TimeIntegration.time(input, state),
                    subdt
                )
            end
        end

        for idx in eachindex(C)
            if C[idx] < zero(Amount)
                C[idx] = zero(Amount)
            end
        end
    end
end


"""
disintegrator(input) -> f!

Prepares the disintegration step. Returns a function `f!(state::State)`. The returned function
modifies the state, popping sediment from the `sediment_buffer` and returns an array of `Amount`.
"""
function disintegrator(input)

    max_h = input.disintegration_rate * input.time.Δt
    w = water_depth(input)
    output = Array{Float64,3}(undef,
                              n_facies(input),
                              input.box.grid_size...)

    depositional_resolution = input.depositional_resolution
    iz = input.intertidal_zone
    tf = input.disintegration_transfer

    return function (state)

        wn = w(state)
        wn .+= iz

        h = min.(max_h, state.sediment_height)
        h[wn .<= 0.0u"m"] .= 0.0u"m"

        state.sediment_height .-= h

        pop_sediment!(state.sediment_buffer,
                      h ./ depositional_resolution .|> NoUnits,
                      output)

        return output .* depositional_resolution
    end
end


"""
transporter(input::Input) -> f

Prepares the transportation step. Returns a function `f(state::State, active_layer)`,
transporting the active layer, returning a transported `Amount` of sediment.
"""
function transporter(input)

    if input.transport_substeps == :adaptive
        return adaptive_transporter(input)
    end

    solver = transport_solver(input.transport_solver, input.box)

    w = water_depth(input)
    box = input.box
    Δt = input.time.Δt / input.transport_substeps
    steps = input.transport_substeps
    fs = input.facies
    iz = input.intertidal_zone

    return function (state)

        wd = w(state)
        wd .+= iz

        C = state.active_layer

        for (i, f) in pairs(fs)
            wv = isnothing(input.wave_velocity) ?
                 f.wave_velocity :
                 input.wave_velocity

            for _ in 1:steps
                solver(
                    (C, _) -> transport(box,
                                        f.diffusion_coefficient,
                                        wv,
                                        C, wd),
                    view(C, i, :, :),
                    TimeIntegration.time(input, state),
                    Δt
                )
            end
        end

        for idx in eachindex(C)
            if C[idx] < zero(Amount)
                C[idx] = zero(Amount)
            end
        end
    end
end


function write_header(input::AbstractInput, output::AbstractOutput)

    set_attribute(output, "intertidal_zone",
                  input.intertidal_zone |> in_units_of(u"m"))
    set_attribute(output, "disintegration_rate",
                  input.disintegration_rate |> in_units_of(u"m/Myr"))

    if input.cementation_time !== nothing
        set_attribute(output, "cementation_time",
                      input.cementation_time |> in_units_of(u"Myr"))
    end
end

end
# ~/~ end