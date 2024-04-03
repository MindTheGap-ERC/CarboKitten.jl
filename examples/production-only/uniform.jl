# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/uniform.jl>>[init]
#| creates: data/ca-prod.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.CaProd
using CarboKitten.BoundaryTrait
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using Unitful

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Periodic{2}}(
    grid_size = (50, 50),
    phys_scale = 1.0u"m"
  ),
  time = TimeProperties(
    Δt = 0.001u"Myr",
    steps = 1000,
    write_interval = 1
  ),
  sea_level=_ -> 0.0u"m",
  subsidence_rate=50.0u"m/Myr",
  initial_depth=_ -> 2.0u"m",
  facies=MODEL1,
  insolation=2000.0u"W/m^2"
)
end

CaProd.main(Script.DEFAULT_INPUT, "data/ca-prod.h5")
# ~/~ end