using Documenter
using CarboKitten

path_to_ext = "C:/Users/lucy.salmon/CarboKitten.jl-main/CarboKitten.jl-main - JossReady/docs_extension"
include(joinpath(path_to_ext, "src", "ExtDocs.jl"))

makedocs(
    sitename = "CarboKitten Extensions",
    modules = [getfield(Main, :ExtDocs)],
    warnonly = [:missing_docs],
    remotes = nothing,
    format = Documenter.HTML(prettyurls = false),
    source = "docs/src",
    build = "docs/build",
    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Workflows" => "workflows.md",
        "Reference" => [
            "Reconstruction" => "reference/reconstruction.md",
            "Visualization" => "reference/visualization.md",
            "Production and compaction" => "reference/production_compaction.md",
            "Environment interpretation" => "reference/environment_interpretation.md",
            "Archive and output" => "reference/archive_output.md",
            "Batch analysis" => "reference/batch_analysis.md",
            "Transport and wave physics" => "reference/transport_wave_physics.md",
        ],
    ]
)