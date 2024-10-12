# ~/~ begin <<docs/src/bosscher-1992.md#examples/bs92/bs92.jl>>[init]
module BS92

using DifferentialEquations
using CSV
using DataFrames
using Interpolations

# ~/~ begin <<docs/src/bosscher-1992.md#b92-model>>[init]
# ~/~ begin <<docs/src/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))
# ~/~ end

struct Parameters
     I₀::Float64
     Iₖ::Float64
     k::Float64
     gₘ::Float64
end

g(p::Parameters, w) = g(p.gₘ, p.I₀, p.Iₖ, p.k, w)
# ~/~ end
# ~/~ begin <<docs/src/bosscher-1992.md#b92-model>>[1]
function model(p::Parameters, s, t_end::Float64, h₀::Float64)
     ∂h(h::Float64, _, t::Float64) = let w = h - s(t)
          w >= 0.0 ? -g(p, h - s(t)) : 0.0
     end
     ode = ODEProblem(∂h, h₀, (0.0, t_end), Nothing)
     solve(ode, Euler(), dt=10.0, reltol=1e-6, saveat=1000.0)
end
# ~/~ end

function sealevel_curve()
     data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
     linear_interpolation(data.time, data.depth)
end

struct Scenario
     param::Parameters
     sealevel
     t_end::Float64
end

model(s::Scenario, h₀::Float64) = model(s.param, s.sealevel, s.t_end, h₀)

SCENARIO_A = Scenario(
     Parameters(2000.0, 250.0, 0.05, 0.005),
     sealevel_curve(),
     80_000.0)

end
# ~/~ end