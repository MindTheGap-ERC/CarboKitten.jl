import CarboKitten.Visualization: resample_to_regular_grid

export resample_to_regular_grid

function resample_to_regular_grid(
    facies_layers::AbstractVector{<:AbstractMatrix{<:Integer}},
    thickness_layers::AbstractVector{<:AbstractMatrix{<:Real}},
    dz::Real
)
    nx, ny = size(facies_layers[1])

    total_thickness = zeros(Float32, nx, ny)
    for k in eachindex(thickness_layers)
        total_thickness .+= Float32.(thickness_layers[k])
    end

    nz = Int(ceil(maximum(total_thickness) / dz))
    grid = zeros(UInt8, nx, ny, nz)

    for i in 1:nx, j in 1:ny
        z_top = 0.0f0

        for k in eachindex(facies_layers)
            h = Float32(thickness_layers[k][i, j])
            f = UInt8(facies_layers[k][i, j])

            z_base = z_top + h

            k1 = Int(floor(z_top / dz)) + 1
            k2 = Int(ceil(z_base / dz))

            for kk in k1:k2
                if kk <= nz
                    grid[i, j, kk] = f
                end
            end

            z_top = z_base
        end
    end

    return grid
end