# ~/~ begin <<docs/src/ca-with-production.md#examples/visualize_production.jl>>[init]
module Script
    using HDF5
    using Plots

    function main()
        h5open("data/ca-prod.h5", "r") do fid
            attr = attributes(fid["input"])
            Δt = attr["delta_t"][]
            subsidence_rate = attr["subsidence_rate"][]
            t_end = fid["input/t"][end-1]
            total_subsidence = subsidence_rate * t_end
            total_sediment = sum(fid["sediment"][:,:,:,:]; dims=3)
            initial_height = fid["input/height"][:]
            elevation = cumsum(total_sediment; dims=4) .* Δt .- initial_height .- total_subsidence
            elevation[1,:,1,:]
        end
    end
end

Script.main()
# ~/~ end