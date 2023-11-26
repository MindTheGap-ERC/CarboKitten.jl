# ~/~ begin <<docs/src/submarine-transport.md#src/CATP.jl>>[init]
module CATP

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.Stencil
using CarboKitten.Burgess2013.CA

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
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-input>>[init]
@kwdef struct Input
  sea_level
  subsidence_rate
  initial_depth

  grid_size::NTuple{2,Int}
  phys_scale::Float64
  Δt::Float64

  time_steps::Int
  write_interval::Int

  facies::Vector{Facies}
  insolation::Float64

  Δz::Float64
  buffer_depth::Int
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
struct TransportFrame
  disintegration::Array{Float64,2}   # x y
  deposition::Array{Float64,3}       # x y f
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[init]
function deposit_material(
    grid_size::NTuple{2, Int},
    Δt::Float64,
    Δz::Float64,
    s::State,
    facies::T) where T <: AbstractArray{Float64,3}

  Threads.@threads for idx in CartesianIndices(grid_size)
    prod = facies[Tuple(idx)..., :]
    Δh = sum(prod) .* Δt
    fractions = prod ./ sum(prod)
    column = s.sediment[Tuple(idx)..., :, :]
    bucket = sum(column[1, :])
    if bucket + Δh / Δz > 1.0
      column[1, :] .+= fractions .* (1.0 - bucket)
      column[2:end, :] = column[1:end-1, :]
      column[1, :] = fractions .* (Δh / Δz - 1.0)
    else
      column[1, :] .+= fractions .* (Δh / Δz)
    end

    s.height[idx] .-= Δh
  end
end

function production_updater(input::Input)
  function (s::State, Δ::ProductFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.production)
  end
end
# ~/~ end
# ~/~ begin <<docs/src/submarine-transport.md#cat-update>>[1]
function transport_updater(input::Input)
  function (s::State, Δ::TransportFrame)
    deposit_material(input.grid_size, input.Δz, s, Δ.deposition)
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