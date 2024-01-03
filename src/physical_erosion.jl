module physical_erosion

using CarboKitten.Burgess2013
using CarboKitten.emperical_denudation

export sediments_redstribution
const kv = 0.23

function D_physical(elevation::Matrix{Float64}, cellsize::Float64, facies::Facies)
    slope = calculate_slope(elevation,cellsize)
    nrows, ncols = size(elevation)
    # function format
    for i in 2:nrows-1
        for j in 2:ncols-1
            D[i,j] = -kv .* (1-facies.inf).^(1/3) .* slope.^(2/3)
        end
    end
    return D
end

function sediments_redistribution(elevation::Matrix{Float64}, cellsize::Float64)
    dh = zeros(Float64,3,3)
    slope = calculate_slope(elevation,cellsize)
    for i in 2:nrows-1
        for j in 2:ncols-1
            (elevation[i - 1, j - 1] .- elevation[i, j]) >= 0 ? dh[i-1,j-1] = D .* ((elevation[i - 1, j - 1] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i-1,j] .- elevation[i, j]) >= 0 ? dh[i-1,j] = sqrt(2) * D .* ((elevation[i-1,j] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i-1,j+1] .- elevation[i, j]) >= 0 ? dh[i-1,j+1] = D .* ((elevation[i-1,j+1] .- elevation[i, j]) ./ cellsize) ./ slope: 0
            (elevation[i,j-1] .- elevation[i, j]) >= 0 ? dh[i,j-1] = sqrt(2) * D .* ((elevation[i,j-1] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i,j+1] .- elevation[i, j]) >= 0 ? dh[i,j+1] = sqrt(2) * D .* ((elevation[i,j+1] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i+1,j-1] .- elevation[i, j]) >= 0 ? dh[i+1,j-1] = D .* ((elevation[i+1,j-1] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i+1,j] .- elevation[i, j]) >=0 ? dh[i+1,j] = sqrt(2) * D .* ((elevation[i+1,j] .- elevation[i, j]) ./ cellsize) ./ slope : 0
            (elevation[i+1,j+1] .- elevation[i, j]) >=0 ? dh[i+1,j+1] = D .* ((elevation[i+1,j+1] .- elevation[i, j]) ./ cellsize) ./ slope : 0
        end
    end
    return dh
end
end