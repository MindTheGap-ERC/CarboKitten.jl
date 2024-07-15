using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "ref.bib"))
# makedocs(; plugins=[bib], ...)

module Entangled
    using DataStructures: DefaultDict

    function transpile_md(src)
        counts = DefaultDict(0)
        Channel{String}() do ch
            for line in src
                if (m = match(r"``` *{[^#}]*#([a-zA-Z0-9\-_]+)[^}]*\}", line)) !== nothing
                    term = counts[m[1]] == 0 ? "≣" : "⊞"
                    put!(ch, "```@raw html")
                    put!(ch, "<div class=\"noweb-label\">⪡" * m[1] * "⪢" * term * "</div>")
                    put!(ch, "```")
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

    function transpile_file(src, target_path)
        mkpath(joinpath(target_path, dirname(src)))
        content = open(readlines, src, "r")
        open(joinpath(target_path, basename(src)), "w") do fout
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

is_markdown(path) = splitext(path)[2] == ".md"
sources = filter(is_markdown, readdir(joinpath(@__DIR__, "src"), join=true))
path = joinpath(@__DIR__, "transpiled")
rm(path; force=true, recursive=true)
mkpath(path)
Entangled.transpile_file.(sources, path)
copydir(joinpath(@__DIR__, "src/fig"), joinpath(path, "fig"))

makedocs(
    source=path,
    sitename="CarboKitten",
    # repo=Remotes.GitHub("MindTheGap-ERC", "CarboKitten"),
    pages = [
        "Introduction" => "index.md",
        "Bosscher and Schlager 1992" => "bosscher-1992.md",
        "CarboCAT" => [
            "Summary" => "carbocat.md",
            "Cellular Automaton" => "carbocat-ca.md",
            "Model with CA and Production" => "ca-with-production.md",
            "Sediment Transport" => "carbocat-transport.md"
        ],
        "Denudation" => [
            "Denudation" => "Denudation.md",
            "Model" => "ca-prod-with-erosion.md"
        ],
        "Algorithms" => [
            "Unitful" => "unitful.md",
            "Boxes" => "boxes.md",
            "Sediment Buffers" => "sediment-buffer.md",
            "Stencils" => "stencils.md",
            "Utility" => "utility.md"
        ],
        "References" => "references.md"
    ],
    plugins = [bib])

deploydocs(
    repo="github.com/MindTheGap-ERC/CarboKitten.jl"
)
