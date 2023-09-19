# ~/~ begin <<docs/src/ca-with-production.md#examples/ca-with-prod.jl>>[init]
using CarboKitten
using CarboKitten.Burgess2013.CA
using CarboKitten.Stencil: Reflected
using CarboKitten.Utility
using CarboKitten.BS92: g

using HDF5

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-input>>[init]
struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

struct Input
    sea_level
    subsidence_rate
    initial_depth

    grid_size::NTuple{2,Int}
    phys_scale::Float64
    Δt::Float64
    time_steps::Int

    facies::Vector{Facies}
    insolation::Float64
end

function production(insolation::Float64, facies::Facies, water_depth::Float64)
    gₘ = facies.maximum_growth_rate
    I₀ = insolation
    Iₖ = facies.saturation_intensity
    w = water_depth
    k = facies.extinction_coefficient
    return gₘ * tanh(I₀/Iₖ * exp(-w * k))
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-frame>>[init]
struct Frame
    production :: Array{Float64, 3}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-state>>[init]
struct State
    time :: Float64
    height :: Array{Float64, 2}
end
# ~/~ end


cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))

function rules(facies::Vector{Facies})
    function (neighbourhood::Matrix{Int}, order::Vector{Int})
        cell_facies = neighbourhood[3, 3]
        neighbour_count(f) = sum(neighbourhood .== f)
        if cell_facies == 0
            for f in order
                n = neighbour_count(species)
                (a, b) = facies[f].activation_range
                if a <= n && n <= b
                    return species
                end
            end
            0
        else
            n = neighbour_count(cell_facies) - 1
            (a, b) = facies[cell_facies].viability_range
            (a <= n && n <= b ? cell_facies : 0)
        end
    end    
end

function run_ca(::Type{B}, facies::Vector{Facies}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    r = rules(facies)
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), r)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[init]
function initial_state(input::Input)  # -> State
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init>>[init]
    # n_facies = length(input.facies)
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(height, 0.0)
    # ~/~ end
end

function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    init = rand(0:n_facies, input.grid_size...)
    ca = run_ca(Reflected{2}, input.facies, init, 3)

    function water_depth(s::State)
        input.sea_level(s.time) - s.height
    end
    # ~/~ end
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(Float64, input.grid_size..., n_facies)
        facies_map = take!(ca)
        w = water_depth(s)
        for idx in CartesianIndices(facies_map)
            if facies_map[idx] == 0
                continue
            end
            result[Tuple(idx)..., f] = production(input.insolation, input.facies[facies_map[idx]], w[idx])
        end
        return result
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[1]
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, d::Frame)
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-update>>[init]
        for f in 1:n_facies
            s.height .+= d[f]
        end
        s.height .-= input.subsidence_rate * input.Δt
        s.time += input.Δt
        # ~/~ end
    end
end
# ~/~ end

function main(input::Input)
    s = initial_state(input)
    p = propagator(input)
    u = updater(input)

    for i in 1:input.time_steps
        Δ = p(s)
        u(s, Δ)
    end
end
# ~/~ end