# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/using_ode.jl>>[init]
module BS92

using CarboKitten.DataSets: bosscher_schlager_1992
using Interpolations
using Unitful

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
     ∂h(h::Float64, t::Float64) = let w = h - s(t)
          w >= 0.0 ? -g(p, h - s(t)) : 0.0
     end

     dt = 1000.0
     times = 0.0:dt:t_end
     result = zeros(Float64, length(times))
     result[1] = h₀
     for (i, t) in enumerate(times[1:end-1])
          h = result[i]
          for j = 0:99
               h += ∂h(h, t + j * dt/100) * (dt/100)
          end
          result[i+1] = h
     end
     return result
end
# ~/~ end

function sealevel_curve()
     data = bosscher_schlager_1992()
     linear_interpolation(data.time / u"yr", data.sealevel / u"m")
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