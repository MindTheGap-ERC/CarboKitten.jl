# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/run.jl>>[init]
#| creates: data/output/bs92.h5

module Script

using CarboKitten
using CarboKitten.DataSets: bosscher_schlager_1992

using Interpolations
using Unitful

function sealevel_curve()
    data = bosscher_schlager_1992()
    linear_interpolation(data.time, data.sealevel)
end

const INPUT = BS92.Input(
    tag = "example model BS92",
    box = Box{Coast}(grid_size=(100, 1), phys_scale=600.0u"m"),
    time = TimeProperties(
      Δt = 10.0u"yr",
      steps = 8000,
      write_interval = 100),
    sea_level = let sc = sealevel_curve()
      t -> -sc(t)
    end,
    initial_topography = (x, y) -> - x / 300.0,
    subsidence_rate = 0.0u"m/yr",
    insolation = 400.0u"W/m^2",
    facies = [BS92.Facies(
      maximum_growth_rate = 0.005u"m/yr",
      saturation_intensity = 50.0u"W/m^2",
      extinction_coefficient = 0.05u"m^-1"
    )])

function main()
    run_model(Model{BS92}, INPUT, "data/output/bs92.h5")
end

end

Script.main()
# ~/~ end
