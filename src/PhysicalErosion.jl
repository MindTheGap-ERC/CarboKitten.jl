module PhysicalErosion

using CarboKitten.Burgess2013
using CarboKitten.Stecil: Periodic
using CarboKitten.EmpericalDenudation
export sediments_redstribution

const kv = 0.23 #very arguable paramster

#calculate how much material could one cell provide to the lower cells
function erosion(cellsize::Float64, Facies::facies)
    slope = slope(cellsize)
    #stencil(Float64,Reflected{2},(3,3),function(w)
    D = -kv .* (1-Facies.inf).^(1/3) .* slope.^(2/3)
    #end)
#=
    nrows, ncols = size(elevation)
    D = similar(elevation)
    # function format
    for i in 2:nrows-1
        for j in 2:ncols-1
            if slope[i,j] <= 0
            D[i,j] = 0
            else
            D[i,j] = -kv .* (1-Facies.inf).^(1/3) .* slope[i,j].^(2/3)
            end
        end
    end
    =#
    return D
end

# find the redistibution co-efficient for the neighboring 8 cells
function sediments_redistribution(cellsize::Float64)
    stencil(Float64,Reflected{2},(3,3),function(w)
    s[1,1] = -(w[1,1] - w[2,2]) / cellsize
    s[1,2] = -(w[1,2] - w[2,2]) / cellsize / sqrt(2)
    s[1,3] = -(w[1,3] - w[2,2]) / cellsize
    s[2,1] = -(w[2,1] - w[2,2]) / cellsize / sqrt(2)
    s[2,2] = -(w[2,2] - w[2,2]) / cellsize
    s[2,3] = -(w[2,3] - w[2,2]) / cellsize / sqrt(2)
    s[3,1] = -(w[3,1] - w[2,2]) / cellsize
    s[3,2] = -(w[3,2] - w[2,2]) / cellsize / sqrt(2)
    s[3,3] = -(w[3,3] - w[2,2]) / cellsize
    sumslope = sum(s)
    s./sumslope
    end)
    #=
    cell_ele = neighbourhood[2,2]
    lower_cells = findall(x -> x < cell_ele, neighbourhood)
    redistribution_matrix = zeros(Float64,3,3)

    if isempty(lower_cells)
        return zeros(Float64, 3, 3)  # No redistribution needed
    end

    slope_matrix = zeros(Float64,3,3)
    slope_matrix[1, 1] = (neighbourhood[1,1] - cell_ele) ./ cellsize
    slope_matrix[1, 2] = (neighbourhood[1,2] - cell_ele) ./ cellsize / sqrt(2)
    slope_matrix[1, 3] = (neighbourhood[1,3] - cell_ele) ./ cellsize
    slope_matrix[2, 1] = (neighbourhood[2,1] - cell_ele) ./ cellsize / sqrt(2)
    slope_matrix[2, 2] = 0
    slope_matrix[2, 3] = (neighbourhood[2,3] - cell_ele) ./ cellsize / sqrt(2)
    slope_matrix[3, 1] = (neighbourhood[3,1] - cell_ele) ./ cellsize
    slope_matrix[3, 2] = (neighbourhood[3,2] - cell_ele) ./ cellsize / sqrt(2)
    slope_matrix[3, 3] = (neighbourhood[3,3] - cell_ele) ./ cellsize

    for i in 1:3
        for j in 1:3
            if slope_matrix[i, j] > 0 
                slope_matrix[i, j] = 0
            else
                slope_matrix[i,j] = -slope_matrix[i,j]
            end
        end
    end
    
    sum_slope = sum(slope_matrix)

    redistribution_matrix .= slope_matrix ./ sum_slope

    return redistribution_matrix
    =#
end

function redistribution_amount(cellsize::Float64,Facies::facies)
    stencil(Float64,Reflected{2},(3,3),function(w)
    s = sediments_redistribution(cellsize)
    D = erosion(cellsize,Facies)
    s.* D[2,2]
    end)
    #=
    nrows, ncols = size(elevation)
    D_matrix = zeros(Float64,3,3)
    result = Matrix{Float64}[]
    for i in 2:nrows-1
        for j in 2:ncols-1
            # Extract the 3x3 neighborhood for the current cell
            m = zeros(Float64,3,3)
            neighbourhood = elevation[i-1:i+1, j-1:j+1]
            m .= sediments_redistribution(neighbourhood,cellsize)
            D = D_physical(neighbourhood, cellsize, Facies)
            D_matrix .= D[2,2] .* m
            # Store the result in the corresponding cell of the results matrix
            push!(result, D_matrix)
        end
    end
    return result
    =#
end

function adderosion(cellsize::Float64, Facies::facies)
    stencil(Float64,Reflected{2},(3,3),function(w)
    amount_matrices = redistribution_amount(cellsize,Facies)
    sum(mount_matrices, dims=(3, 4))
    end)
end

end