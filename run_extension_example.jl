import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "examples", "extension", "run_extension_example.jl"))