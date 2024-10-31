# ~/~ begin <<docs/src/cases/tabular-sea-level.md#examples/tabular-sea-level.jl>>[init]
module TabularSeaLevel

using DelimitedFiles: readddlm
using DataFrames
using CarboKitten.DataSets: artifact_dir

# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[init]
using CarboKitten.Components.Common
using CarboKitten.Model.ALCAP2
# ~/~ end
# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[1]
function miller_2020()
    dir = artifact_dir()
    filename = joinpath(dir, "Miller2020", "Cenozoic_sea_level_reconstruction.tab")

    data, header = readdlm(filename, '\t', header=true)
    return DataFrame(
        time=-data[:,4] * u"kyr",
        sealevel=data[:,7] * u"m",
        reference=categorical(data[:,3]))
end
# ~/~ end

end

TabularSeaLevel.main()
# ~/~ end
