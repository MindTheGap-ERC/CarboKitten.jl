using Documenter
using DocumenterCitations
using DocumenterMermaid

using CarboKitten

bib = CitationBibliography(joinpath(@__DIR__, "src", "ref.bib"))
# makedocs(; plugins=[bib], ...)

include("component_graphs.jl")

module Entangled
using DataStructures: DefaultDict

function transpile_md(src)
    counts = DefaultDict(0)
    Channel{String}() do ch
        for line in src
            if (m = match(r"( *)``` *{[^#}]*#([a-zA-Z0-9\-_]+)[^}]*\}", line)) !== nothing
                term = counts[m[2]] == 0 ? "≣" : "⊞"
                put!(ch, "$(m[1])```@raw html")
                put!(ch, "$(m[1])<div class=\"noweb-label\">⪡" * m[2] * "⪢" * term * "</div>")
                put!(ch, "$(m[1])```")
                put!(ch, line)
                counts[m[1]] += 1
            elseif (m = match(r"``` *{[^}]*file=([a-zA-Z0-9\-_\.\/\\]+)[^}]*}", line)) !== nothing
                put!(ch, "```@raw html")
                put!(ch, "<div class=\"noweb-label\">file:<i>" * m[1] * "</i></div>")
                put!(ch, "```")
                put!(ch, line)
            else
                put!(ch, line)
            end
        end
    end
end

function transpile_file(srcdir, file, target_path)
    mkpath(joinpath(target_path, dirname(file)))
    content = open(readlines, joinpath(srcdir, file), "r")
    open(joinpath(target_path, file), "w") do fout
        join(fout, transpile_md(content), "\n")
    end
end
end

function copydir(src, dst)
    for (path, subdir, files) in walkdir(src)
        mkpath(joinpath(dst, subdir...))
        for f in files
            cp(joinpath(path, subdir..., f), joinpath(dst, subdir..., f))
        end
    end
end

function get_source_files()
    home = joinpath(@__DIR__, "src")
    is_markdown(path) = splitext(path)[2] == ".md"
    files_in_dir(d, _, fs) = (relpath(joinpath(d, f), home) for f in filter(is_markdown, fs))
    Iterators.flatten(Iterators.map(
        splat(files_in_dir), walkdir(home))) |> collect
end

sources = get_source_files()
path = joinpath(@__DIR__, "transpiled")
rm(path; force=true, recursive=true)
mkpath(path)
Entangled.transpile_file.(joinpath(@__DIR__, "src"), sources, path)
run(`touch $(joinpath(path, "first_tutorial.md"))`)
copydir(joinpath(@__DIR__, "src/fig"), joinpath(path, "fig"))

makedocs(
    source=path,
    modules = [CarboKitten],
    sitename="CarboKitten",
    # repo=Remotes.GitHub("MindTheGap-ERC", "CarboKitten"),
    pages=[
        "Introduction" => "index.md",
        "Models" => [
            "Bosscher and Schlager 1992" => "bosscher-1992.md",
            "Model with CA and Production" => "ca-with-production.md",
            # "With Denudation" => "ca-prod-with-denudation.md",
            "ALCAPS" => "model-alcap.md"
        ],
        "Examples" => [
            "Tutorial (Pluto notebook)" => "first_tutorial.md",
            "Tabular Sea Levels" => "cases/tabular-sea-level.md"
        ],
        "Architecture" => "architecture.md",
        "Model Components" => [
            "Components" => "components/components.md",
            "Tags" => "components/tag.md",
            "Boxes" => "components/boxes.md",
            "Time" => "components/time.md",
            "Facies" => "components/facies.md",
            "Cellular Automata" => "components/cellular-automata.md",
            "Water Depth" => "components/waterdepth.md",
            "Production" => "components/production.md",
            "HDF5 Writer" => "components/hdf5.md",
            "Sediment Buffers" => "components/sediment_buffer.md",
            "Active Layer Transport" => "active-layer-transport.md",
            "Onshore Transport" => "onshore-transport.md"
        ],
        "Visualizations" => "visualization.md",
        "CarboCAT" => [
            "Summary" => "carbocat.md",
        ],
        "Denudation" => [
            "Denudation" => "denudation/denudation.md",
            "Empirical Denudation" => "denudation/empirical.md",
            "Chemical Dissolution" => "denudation/chemical.md",
            "Physical Erosion" => "denudation/physical_erosion.md"
        ],
        "Input & Output" => [
            "Input Methods" => "input-methods.md",
            "CSV Export" => "data-export.md"
        ],
        "Algorithms" => [
            "Unitful" => "unitful.md",
            "Boxes" => "boxes.md",
            "Stencils" => "stencils.md",
            "Utility" => "utility.md"
        ],
        "API Documentation" => "api.md",
        "References" => "references.md"
    ],
    plugins=[bib])

cp(joinpath(@__DIR__, "notebooks/first_tutorial.html"), joinpath(@__DIR__, "build/first_tutorial/index.html"), force=true)

deploydocs(
    repo="github.com/MindTheGap-ERC/CarboKitten.jl"
)
