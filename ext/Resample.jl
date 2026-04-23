import CarboKitten.Visualization:resample_to_regular_grid
     

export resample_to_regular_grid

function resample_to_regular_grid(
    facies_layers::Vector{Matrix{UInt8}},
    thickness_layers::Vector{Matrix{Float32}},
    dz::Float32
)
    nx, ny = size(facies_layers[1])

    # total thickness per column
    total_thickness = zeros(Float32, nx, ny)
    for k in eachindex(thickness_layers)
        total_thickness .+= thickness_layers[k]
    end

    nz = Int(ceil(maximum(total_thickness) / dz))
    grid = zeros(UInt8, nx, ny, nz)

    # loop over columns
    for i in 1:nx, j in 1:ny

        z_top = 0.0f0

        for k in eachindex(facies_layers)
            h = thickness_layers[k][i,j]
            f = facies_layers[k][i,j]

            z_base = z_top + h

            # find voxel range intersecting this layer
            k1 = Int(floor(z_top / dz)) + 1
            k2 = Int(ceil(z_base / dz))

            for kk in k1:k2
                if kk <= nz
                    grid[i,j,kk] = f
                end
            end

            z_top = z_base
        end
    end

    return grid
end