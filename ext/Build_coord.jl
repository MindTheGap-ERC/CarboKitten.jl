import CarboKitten.Visualization:
     build_vertical_coordinates

export build_vertical_coordinates

function build_vertical_coordinates(thickness_layers)
    # thickness_layers: Vector of (nx, ny) matrices

    nx, ny = size(thickness_layers[1])
    nl = length(thickness_layers)

    z_top = [zeros(Float32, nx, ny) for _ in 1:nl]
    z_base = [zeros(Float32, nx, ny) for _ in 1:nl]

    for i in 1:nx, j in 1:ny
        z = 0.0f0
        for k in reverse(1:nl)  # bottom → top
            h = thickness_layers[k][i,j]

            z_top[k][i,j] = z
            z_base[k][i,j] = z + h

            z += h
        end
    end

    return z_top, z_base
end