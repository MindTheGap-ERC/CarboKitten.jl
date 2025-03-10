module TestexportSpec

using Test
using DataFrames
using CSV

const PATH = "data/output"
const TAG = "alcap-example-deposition_only"

data = CSV.read("$(PATH)/$(TAG)_sc.csv", DataFrame)
sc1_f1 = data.var"sc1_f1 [m]"
sc1_f2 = data.var"sc1_f2 [m]"
sc1_f3 = data.var"sc1_f3 [m]"
sum(sc1_f1+sc1_f2+sc1_f3)

data2 = CSV.read("$(PATH)/$(TAG)_sac.csv", DataFrame)
sac1 = data2.var"sac1 [m]"[end]

@test sum(sc1_f1+sc1_f2+sc1_f3) ≈ sac1
end