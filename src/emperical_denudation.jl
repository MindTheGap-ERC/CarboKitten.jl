## this code used the emperical equations with form of D = f(precipitation, slope) to estimate the denudation rates on exposed carbonate platform.
module emperical_denudation

using CarboKitten.CaProd
using CarboKitten.Stecil: Periodic
export empericaldenudation
# calculate planar slopes based on [ARCGIS apporach](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm)


function slope(cellsize::Float64)
    stencil(Float64,Reflected{2},(3,3),function(w)
    dzdx = (-w[1,1] - 2 * w[2,1] -w[3,1] + w[1,3] + 2 * w[2,3] + w[3,3])/(8*cellsize)
    dzdy = (-w[1,1] - 2 * w[1,2] -w[1,3] + w[3,1] + 2 * w[3,2] + w[1,1])/(8*cellsize)
    atan(sqrt(dzdx^2 + dzd^2))  * (180 / π)
    end)
end



function empericaldenudation(precip::Float64, water_depth::Matrix{Float64}, cellsize::Float64)
    const a = 9.1363
    const b = -0.008519
    const c = 580.51
    const d = 9.0156
    const e = -0.1245
    const f = 4.91086
    if water_depth < 0
    slope = slope(cellsize)
    D .= (a ./ (1 .+ exp.(b.*(precip .- c)))) .* (d ./ (1 .+ exp.(e.*(slope .- f)))) # using logistic function
    end
    return D./1000 #m/kyr
end

end