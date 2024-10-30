# ~/~ begin <<docs/src/utility.md#src/DataSets.jl>>[init]
module DataSets

export read_tsv, read_csv

using DelimitedFiles: readdlm
using DataFrames
using CategoricalArrays
using CSV
using Pkg.Artifacts
using Unitful

function read_tsv(filename)
    data, header = readdlm(filename, '\t', header=true)
    return DataFrame(data, vec(header))
end

function read_csv(filename)
    return DataFrame(CSV.File(filename))
end

function artifact_dir()
    subfolder = first(readdir(artifact"data"))
    return joinpath(artifact"data", subfolder)
end

function bosscher_schlager_1992()
    dir = artifact_dir()
    filename = joinpath(dir, "Bosscher1992", "bs92-sealevel-curve.csv")
    df = read_csv(filename)
    return DataFrame(time=df.time * u"yr", sealevel=df.depth * u"m")
end

function miller_2020()
    dir = artifact_dir()
    filename = joinpath(dir, "Miller2020", "Cenozoic_sea_level_reconstruction.tab")
    df = read_tsv(filename)
    return DataFrame(
        time=-df[!,4] * u"kyr",
        sealevel=df[!,7] * u"m",
        reference=categorical(df[!,2]))
end

end
# ~/~ end
