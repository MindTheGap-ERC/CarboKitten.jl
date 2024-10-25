# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/multi-facies-run.jl>>[init]
#| creates: data/output/bs92-multi-facies.h5
#| requires: data/bs92-sealevel-curve.csv

module Script

using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Model.BS92

using Unitful

const FACIES = [
    BS92.Facies(
         maximum_growth_rate=500u"m/Myr"/4,
         extinction_coefficient=0.8u"m^-1",
         saturation_intensity=60u"W/m^2"),
    BS92.Facies(
         maximum_growth_rate=400u"m/Myr"/4,
         extinction_coefficient=0.1u"m^-1",
         saturation_intensity=60u"W/m^2"),
    BS92.Facies(
         maximum_growth_rate=100u"m/Myr"/4,
         extinction_coefficient=0.005u"m^-1",
         saturation_intensity=60u"W/m^2")]
	
const INPUT = BS92.Input(
    tag = "example model BS92",
    box = Common.Box{Shelf}(grid_size=(100, 1), phys_scale=150.0u"m"),
    time = TimeProperties(
        Δt = 200.0u"yr",
        steps = 5000,
        write_interval = 1),
    sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    bedrock_elevation = (x, y) -> - x / 300.0,
    subsidence_rate = 50.0u"m/Myr",
    insolation = 400.0u"W/m^2",
    facies = FACIES)

function main()
    H5Writer.run(Model{BS92}, INPUT, "data/output/bs92-multi-facies.h5")
end

end

Script.main()
# ~/~ end
