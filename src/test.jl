#=
###test chemical dissolution###
const T = 288
const P = 1000
const inf = 0.5
const density = 2730

const A = -0.4883+8.074*0.0001*(T-273)
const B = -0.3241+1.6*0.0001*(T-273)
const IA = 0.1 # ion activity
const pco2 = 10^(-1.5)
# chemical basic parameters
const K1::Float64 = 10^(-356.3094-0.06091964*T+21834.37/T+126.8339*log10(T)-1684915/(T^2))
const K2::Float64 = 10^(-107.881-0.03252849*T+5151.79/T+38.92561*log10(T)-563713.9/(T^2))
const KC::Float64 = 10^(-171.9065-0.077993*T+2839.319/T+71.595*log10(T))
const KH::Float64 = 10^(108.3865+0.01985076*T-6919.53/T-40.4515*log10(T)+669365/(T^2))
const gama_Ca::Float64 = 10^(-4A*sqrt(IA)/(1+10^(-8)*B*sqrt(IA)))
const gama_alk::Float64 = 10^(-A*sqrt(IA)/(1+5.4*10^(-8)*B*sqrt(IA)))

# output

ceq = ((pco2) .* (K1 .* KC .* KH) ./(4 * K2 .* gama_Ca .* (gama_alk).^2)).^(1/3)
Deq =  P .* inf * 40 * 1000 * ceq ./ density  # in mm/yr or m/kyr

println(ceq,Deq)
const precip = 1000
const inf = 0.5
const alpha = 2e-6
const L = 10
const density = 2730
const z0 = 10
lambda = precip .* inf ./ (alpha .* L)
D = (precip .* inf .* ceq ./density) .* (1 - (lambda./z0).* (1 - exp(-z0./lambda)))
print(D)=#
using CarboKitten.Stencil

function calculate_slope(elevation::Matrix{Float64}, cellsize::Float64) 
    nrows, ncols = size(elevation)
    slope = similar(elevation)
    #slope_dir = zeros(Float64,size(elevation))
    for i in 2:nrows-1
        for j in 2:ncols-1
            dzdx = ((elevation[i - 1, j + 1] + 2*elevation[i, j + 1] + elevation[i + 1, j + 1]) - (elevation[i - 1, j - 1] + 2*elevation[i, j - 1] + elevation[i + 1, j - 1])) ./ (8 * cellsize)
            dzdy = ((elevation[i + 1, j - 1] + 2*elevation[i + 1, j] + elevation[i + 1, j + 1]) - (elevation[i - 1, j - 1] + 2*elevation[i - 1, j] + elevation[i - 1, j + 1]))/ (8 * cellsize)
            slope[i, j] = atan.(sqrt.(dzdx.^2 + dzdy.^2))  * (180 / π)
            #slope_dir[i,j] = atan.(dzdy,dzdx)
        end
    end

    return slope
    #return slope_dir
end



function slope(cellsize::Float64)
    stencil(Float64,Reflected{2},(3,3),function(w)
    dzdx = (-w[1,1] - 2 * w[2,1] -w[3,1] + w[1,3] + 2 * w[2,3] + w[3,3])/(8*cellsize)
    dzdy = (-w[1,1] - 2 * w[1,2] -w[1,3] + w[3,1] + 2 * w[3,2] + w[1,1])/(8*cellsize)
    atan(sqrt(dzdx^2 + dzd^2))  * (180 / π)
    end)
end

function slope2(cellsize::Float64)
    kernel = [-1 0 1; -2 0 2; -1 0 1] ./ (8*cellsize)
    stencil(Float64, Reflected{2}, (3,3), function (w)
        dzdx = sum(w .* kernel)
        dzdy = sum(w .* kernel')
        atan(sqrt(dzdx^2 + dzd^2))  * (180 / π)
    end)
end

#=
function calculate_D(precip::Float64, elevation::Matrix{Float64}, cellsize::Float64)
    slope = calculate_slope(elevation,cellsize)
    nrows, ncols = size(elevation)

    D = similar(elevation)
    # function format
    for i in 2:nrows-1
        for j in 2:ncols-1
    D[i,j] = (9.1363 ./ (1 .+ exp.(-0.008519.*(precip .- 580.51)))) .* (9.0156 ./ (1 .+ exp.(-0.1245.*(slope[i,j] .- 4.91086))))
        end
    end
    return D
end
=#



const kv = 0.23 #very arguable paramster

#calculate how much material could one cell provide to the lower cells
function D_physical(elevation::Matrix{Float64}, cellsize::Float64, inf::Float64)
    slope = calculate_slope(elevation,cellsize)
    nrows, ncols = size(elevation)
    D = similar(elevation)
    # function format
    for i in 2:nrows-1
        for j in 2:ncols-1
            if slope[i,j] <= 0
            D[i,j] = 0
            else
            D[i,j] = -kv .* (1-inf).^(1/3) .* slope[i,j].^(2/3)
            end
        end
    end
    return D
end

# find the redistibution co-efficient for the neighboring 8 cells
function sediments_redistribution(neighbourhood::Matrix{Float64},cellsize::Float64)
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

end

function calculate_redistribution(elevation::Matrix{Float64},cellsize::Float64,inf::Float64)
    nrows, ncols = size(elevation)
    D_matrix = zeros(Float64,3,3)
    result = Matrix{Float64}[]
    for i in 2:nrows-1
        for j in 2:ncols-1
            # Extract the 3x3 neighborhood for the current cell
            m = zeros(Float64,3,3)
            neighbourhood = elevation[i-1:i+1, j-1:j+1]
            m .= sediments_redistribution(neighbourhood,cellsize)
            D = D_physical(neighbourhood, cellsize, inf)
            D_matrix .= D[2,2] .* m
            # Store the result in the corresponding cell of the results matrix
            push!(result, D_matrix)
        end
    end
    return result
end


#function 


matrix = rand(Float64, 5, 5) * 3
matrix[2, 2] = 1

ele = matrix
csz = 10.0
slope = calculate_slope(ele,csz)
slope_dir = calculate_slope(ele,csz)
D_result = D_physical(ele, csz, 0.5)
RS = sediments_redistribution(ele,csz)
a = calculate_redistribution(ele,csz,0.5)
println(a)