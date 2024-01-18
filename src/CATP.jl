# ~/~ begin <<docs/src/submarine-transport.md#src/CATP.jl>>[init]
module CATP

using ..Stencil: Periodic
using ..Utility
using ..Stencil
using ..Burgess2013.CA
using ..SedimentStack

using HDF5

using .Iterators: drop, peel, partition, map, take
# ~/~ begin <<docs/src/submarine-transport.md#cat-facies>>[init]
struct Facies
  viability_range::Tuple{Int,Int}
  activation_range::Tuple{Int,Int}

  maximum_growth_rate::Float64
  extinction_coefficient::Float64
  saturation_intensity::Float64

  grain_size::Float64
  excess_density::Float64
  critical_angle::Float64
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-input>>[init]
@kwdef struct Input
  sea_level
  subsidence_rate
  initial_depth

  grid_size::NTuple{2,Int}
  boundary::Boundary{2}
  phys_scale::Float64     # km / pixel
  Δt::Float64             # Ma / timestep ?

  time_steps::Int
  write_interval::Int

  facies::Vector{Facies}
  insolation::Float64

  Δz::Float64
  buffer_depth::Int
  # To estimate the average grain size and overdensity of the sediment (and
  # thereby the shear stress), we need to peek at the composition of the top
  # layer, the `top_layer_thickness` parameter controls how deep we sample.
  top_layer_thickness::Float64
  # Wave shear stress function, is a function of depth.
  wave_shear_stress
  # Gravitational accelleration
  g::Float64
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-state>>[init]
mutable struct State
  time::Float64
  height::Array{Float64,2}            # x y
  sediment::Array{Float64,4}          # x y z (f..., t)
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-frame>>[init]
struct ProductFrame
  production::Array{Float64,3}        # x y f
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-frame>>[1]
struct RemoveFrame
  disintegration::Array{Float64,2}   # x y
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[init]
function deposit_material(
    grid_size::NTuple{2, Int},
    Δz::Float64,
    s::State,
    material::AbstractArray{Float64,3})

  Threads.@threads for idx in CartesianIndices(grid_size)
    production = material[Tuple(idx)..., :] .* (1.0 / Δz)
    push_sediment!(s.sediment, production)
    s.height[idx] .-= sum(material[Tuple(idx)..., :]) * Δt
  end
end

function remove_material(
    grid_size::NTuple{2, Int},
    n_facies::Int,
    Δz::Float64,
    s::State,
    thickness::AbstractArray{Float64,2})

  material = Array{Float64,3}(undef, grid_size..., n_facies)
  Threads.@threads for idx in CartesianIndices(grid_size)
    material[Tuple(idx)..., :] = pop_sediment!(s.sediment, thickness[idx])
  end
  return material
end

function deposit_updater(input::Input)
  function (s::State, Δ::ProductFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.production)
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[1]
function remove_updater(input::Input)
  function (s::State, Δ::RemoveFrame)
    loose = remove_material(input.grid_size, length(input.facies), input.Δz, s, Δ.disintegration)
    # need physics to disperse loose material
    return ProductFrame(loose)
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[2]
function time_updater(input::Input)
  function (s::State, Δ::Nothing = nothing)
    s.height .+= input.subsidence_rate * input.Δt
    s.time += input.Δt
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-propagator>>[init]
Vec2 = @NamedTuple{x::Float64, y::Float64}
Base.abs2(a::Vec2) = a.x^2 + a.y^2
Base.abs(a::Vec2) = √(abs2(a))

function shear_stress(input::Input)
  # To compute critical shear stress, take weighted average of critical angles
  θ_f = [sin(f.critical_angle) for f in input.facies]
  gradient_stencil = stencil(Float64, Vec2, input.boundary, (3, 3), function (w)
    kernel = [-1 0 1; -2 0 2; -1 0 1] ./ (8 * phys_scale)
    ( x = sum(kernel .* w)
    , y = sum(kernel' .* w) )
  end)

  ∇ = Matrix{Vec2}(undef, input.grid_size...)
  ρD_f = [f.excess_density * f.grain_size for f in input.facies]

  function stress(s::State)
    gradient_stencil(s.height, ∇)
    α = atan.(abs.(∇))

    top_sediment = peek_sediment(s.sediment, input.top_layer_thickness)
    A = sum(top_sediment .* ρD_f .* input.g; dims=3) ./ sum(top_sediment)
    θ = sum(top_sediment .* θ_f; dims=3) ./ sum(top_sediment) 

    gravitational_stress = A .* sin.(α)
    wave_stress = input.wave_shear_stress.(s.height)
    critical_stress = A .* θ
  end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagator>>[init]
function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # ~/~ end
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(Float64, input.grid_size..., n_facies)
        facies_map, ca = peel(ca)
        w = water_depth(s)
        Threads.@threads for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
            if f == 0
                continue
            end
            result[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
        end
        return Frame(result)
        # ~/~ end
    end
end
# ~/~ end

end  # CATP
# ~/~ end