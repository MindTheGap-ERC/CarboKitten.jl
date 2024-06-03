# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/caps-osc.jl>>[init]
#| creates: data/caps-osc.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using CarboKitten.CaProd
using Unitful

const PERIOD = 200.0u"kyr"
const AMPLITUDE = 4.0u"m"

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Shelf}(
    grid_size = (100, 50),
    phys_scale = 0.15u"km"
    # equivalent:
    # phys_scale = 150"m"
  ),
  time = TimeProperties(
    Δt = 0.0001u"Myr",
    # equivalent: 
    # Δt = 1u"kyr",
    steps = 10000,
    write_interval = 10
  ),
  sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD), 
  subsidence_rate=50.0u"m/Myr",
  initial_depth=x -> x / 300.0,
  facies=MODEL1,
  insolation=400.0u"W/m^2"
)
end

Script.CaProd.main(Script.DEFAULT_INPUT, "data/caps-osc.h5")
# ~/~ end