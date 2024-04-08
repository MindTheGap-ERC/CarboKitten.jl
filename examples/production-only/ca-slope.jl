# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/ca-slope.jl>>[init]
#| creates: data/ca-prod-slope.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using CarboKitten.CaProd
using Unitful

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Shelf}(
    grid_size = (100, 50),
    phys_scale = 1.0u"km"
  ),
  time = TimeProperties(
    Δt = 0.001u"Myr",
    steps = 1000,
    write_interval = 1
  ),
  sea_level=_ -> 0.0u"m",
  subsidence_rate=50.0u"m/Myr",
  initial_depth=x -> x / 2000.0,
  facies=MODEL1,
  insolation=2000.0u"W/m^2"
)
end

Script.CaProd.main(Script.DEFAULT_INPUT, "data/ca-prod-slope.h5")
# ~/~ end