# ~/~ begin <<docs/src/submarine-transport.md#src/CATP.jl>>[init]
module CATP

using ..Vectors
using ..Stencil: Periodic
using ..Utility
using ..Stencil
using ..Burgess2013.CA
using ..SedimentStack
using ..Transport

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
  # critical_angle::Float64

  # magnitude of critical stress, should be ΔρDg sin α where α is the critical angle
  critical_stress::Float64
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
  # wave shear stress function, is a function of depth, and should return a Vec2
  wave_shear_stress
  # Gravitational accelleration
  g::Float64
  # for the transport step, parcels of sediment are converted into particles
  # by setting this parameter to something >1 you can control how many particles
  # are created on each axis. So number of particles is grid width * height * transport_subsample^2.
  transport_subsample::Int
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
  production::Array{Float64,3}        # facies x y
end

# struct Properties
#   excess_density::Float64
#   grain_size::Float64
# end

const Particle = Transport.Particle{Nothing}
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-frame>>[1]
struct ReductionFrame
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
  function (s::State, Δ::ReductionFrame)
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
function submarine_transport(input::Input)
  phys_size = (x = input.grid_size[1] * input.phys_scale,
               y = input.grid_size[2] * input.phys_scale)
  box = Transport.Box(input.grid_size, phys_size, input.phys_scale, length(input.facies))

  # This function does not modify the state, rather it transports the
  # sediments in a given product frame and gives a new product frame
  # as output.
  function (s::State, Δ::ProductFrame)  # -> ProductFrame
    output = Array{Float64}(undef, size(Δ.production))
    transport = Transport.transport(input.boundary, box, function (p::Particle)
      z, ∇ = Transport.interpolate(input.boundary, box, s.height, p.position)
      α = atan(abs(∇))
      Ĝ = ∇ / abs(∇)

      τ_wave = input.wave_shear_stress(z)

      Δρ = input.facies[p.facies].excess_density
      D = input.facies[p.facies].grain_size
      g = input.g
      τ_grav = (Δρ * D * g * sin(α)) * Ĝ

      return τ_grav + τ_wave
    end)
    deposit = Transport.deposit(input.boundary, box, output)

    # iterate all particles, transport them and plot
    for (i, (idx, mass)) in enumerate(pairs(Δ.production))
      subgrid_spacing = 1.0 / input.transport_subsample
      subgrid_axis = 0.0:subgrid_spacing:1.0 - subgrid_spacing
      subgrid = product(subgrid_axis, subgrid_axis)
      for (j, dx) in enumerate(subgrid)
        facies_type = idx[1]
        p = (x=(idx[2] + dx[1]) * input.phys_scale, y=(idx[3] + dx[2]) * input.phys_scale)
        θ = input.facies[facies_type].critical_stress
        particle = Particle(p, mass * subgrid_spacing^2, θ, facies_type, nothing) |> transport |> deposit
      end
    end
    return ProductFrame(output)
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