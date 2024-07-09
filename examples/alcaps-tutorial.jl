# # CarboKitten tutorial

# CarboKitten uses [`Unitful`](@ref unitful) quantities throughout.

using Unitful

# All our plots are made with [`Makie.jl`](https://docs.makie.org/stable/).

using CairoMakie

using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Model.ALCAPS: Input, FACIES, run_model
using CarboKitten.Visualization: sediment_profile!

# ## Config type
#
# The `Config` type is like a dynamic object. Not sure if this is the way to go.
# If we proceed with this, it should go into `CarboKitten.Utility`.

struct Config
	dict::IdDict
	Config(;kwargs...) = new(IdDict(kwargs...))
end

Base.getproperty(cfg::Config, sym::Symbol) = getfield(cfg, :dict)[sym]

# A large part of the input config is the configuration per facies. We'll skip that for the moment and
# start with a default set of facies.

# ## Running a complete example

const INPUT = Input(
    box = Box{Shelf}(grid_size=(50,100), phys_scale=100u"m"),
    time = TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1),

    ca_interval = 1,

    bedrock_elevation        = (x, y) -> -x / 300.0 ,
    sea_level                = t -> 0.0u"m",
    subsidence_rate          = 50.0u"m/Myr",
    disintegration_rate      = 500.0u"m/Myr",
    insolation               = 400.0u"W/m^2",
    sediment_buffer_size     = 50,
    depositional_resolution  = 0.5m,
    facies                   = FACIES))

run_model(INPUT, "tutorial1.h5")

# Next we may visualize the output.

let
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, "tutorial1.h5")
    fig
end
