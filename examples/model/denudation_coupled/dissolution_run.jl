module Script

using Unitful
using CarboKitten
using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Components.Denudation
using CarboKitten.Models: WithDenudation as WDn
using CarboKitten.Export: data_export, CSV, read_slice
using CarboKitten: Box, TimeProperties, OutputSpec, run_model, Model
using CarboKitten.Denudation: Dissolution, EmpiricalDenudation, PhysicalErosion, NoDenudation
using .DenudationParamConfig

const PATH = "data/output"
const TAG = "denudation"
const SURF = [10u"m^2/m^3", 10u"m^2/m^3", 10u"m^2/m^3"]
const DENS = [2730u"kg/m^3", 2730u"kg/m^3", 2730u"kg/m^3"]
const INF = [0.5, 0.5, 0.5]
const DENUDATION = Dissolution(temp=293.0u"K", precip=1.0u"m", pco2=10^(-2.5) * u"atm", reactionrate=2e-3u"m/yr")

function main()
    FACIES = facies(SURF, DENS, INF)
    INPUT = input(TAG, DENUDATION, FACIES)
    run_model(Model{WDn}, INPUT, "$(PATH)/$(TAG).h5")
    header, profile = read_slice("$(PATH)/$(TAG).h5", :profile)
    columns = [profile[i] for i in 10:20:70]
    data_export(
        CSV(:sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
            :water_depth => "$(PATH)/$(TAG)_wd.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
         header,
         columns)
end
end

Script.main()
