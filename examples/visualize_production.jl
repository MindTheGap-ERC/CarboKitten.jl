# ~/~ begin <<docs/src/ca-with-production.md#examples/visualize_production.jl>>[init]
module Script
    using HDF5
    using Plots

    function main()
        h5open("data/ca-prod.h5", "r") do fid
            total_sediment = sum(fid["sediment"][:,:,:,:]; dims=3)
            initial_height = fid["input"]["height"][:]
            height = initial_height' .- cumsum(total_sediment; dims=4)
            height[1,:,1,:]
        end
    end
end

Script.main()
# ~/~ end