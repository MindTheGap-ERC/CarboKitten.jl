## this code used the emperical equations with form of D = f(precipitation, slope) to estimate the denudation rates on exposed carbonate platform.
#module EmpericalDenudation

using CarboKitten.Stencil: Periodic, stencil
export emperical_denudation,slope_kernel

#calculate planar slopes based on [ARCGIS apporach](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm)

function slope_kernel(w::Matrix{Float64},cellsize::Float64)
    dzdx = (-w[1,1] - 2 * w[2,1] -w[3,1] + w[1,3] + 2 * w[2,3] + w[3,3])/(8*cellsize)
    dzdy = (-w[1,1] - 2 * w[1,2] -w[1,3] + w[3,1] + 2 * w[3,2] + w[1,1])/(8*cellsize)
    atan(sqrt(dzdx^2 + dzdy^2))  * (180 / π)
end

slopefn(w,slope_matrix,cellsize::Float64)

const a = 9.1363
const b = -0.008519
const c = 580.51
const d = 9.0156
const e = -0.1245
const f = 4.91086

# calculate denudation based on regressed function
function emperical_denudation(precip::Float64, slope::Float64)
    (a ./ (1 .+ exp.(b.*(precip .- c)))) .* (d ./ (1 .+ exp.(e.*(slope .- f)))) # using logistic function
end


#=
function emperical_denudation(precip::Float64, cellsize::Float64)

    stencil_result = stencil(Float64, Periodic{2}, (3, 3), (w) -> begin
        dzdx = (-w[1, 1] - 2 * w[2, 1] - w[3, 1] + w[1, 3] + 2 * w[2, 3] + w[3, 3]) / (8 * cellsize)
        dzdy = (-w[1, 1] - 2 * w[1, 2] - w[1, 3] + w[3, 1] + 2 * w[3, 2] + w[1, 1]) / (8 * cellsize)
        atan(sqrt(dzdx^2 + dzdy^2)) * (180 / π)
    end)
    return stencil_result
    #(a ./ (1 .+ exp.(b .* (precip .- c)))) .* (d ./ (1 .+ exp.(e .* (stencil_result .- f))))
end

function emperical_denuadtion(precip::Float64,cellsize::Float64)
    
    (a ./ (1 .+ exp.(b .* (precip .- c)))) .* (d ./ (1 .+ exp.(e .* (stencil_result .- f))))
end
=#
w = rand(50,100)
result = emperical_denudation(1000.0,1.0)

#end