# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/run.jl>>[init]
#| creates: data/output/bs92.h5
#| requires: data/bs92-sealevel-curve.csv

module Script
# using Logging
# using TerminalLoggers
# global_logger(TerminalLogger(right_justify=80))

using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Model.BS92

using CSV
using DataFrames
using Interpolations
using Unitful

function sealevel_curve()
    data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
    linear_interpolation(data.time, data.depth)
end

const INPUT = Input(
    tag = "example model BS92",
    box = Common.Box{Shelf}(grid_size=(100, 1), phys_scale=600.0u"m"),
    time = TimeProperties(
      Δt = 10.0u"yr",
      steps = 8000,
      write_interval = 100),
    sea_level = let sc = sealevel_curve()
      t -> -sc(t / u"yr") * u"m"
    end,
    bedrock_elevation = (x, y) -> - x / 300.0,
    subsidence_rate = 0.0u"m/yr",
    insolation = 400.0u"W/m^2",
    facies = [Facies(
      maximum_growth_rate = 0.005u"m/yr",
      saturation_intensity = 50.0u"W/m^2",
      extinction_coefficient = 0.05u"m^-1"
    )])

function main()
    H5Writer.run(Model{BS92}, INPUT, "data/output/bs92.h5")
end

end

Script.main()
# ~/~ end