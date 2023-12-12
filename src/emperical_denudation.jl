## this code used the emperical equations with form of D = f(precipitation, slope) to estimate the denudation rates on exposed carbonate platform.
module emperical_denudation

using CarboKitten.CaProd
export D
# calculate planar slopes based on [ARCGIS apporach](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm)


function calculate_slope(elevation::Matrix{T}, cellsize::T) where T
    nrows, ncols = size(elevation)
    slope = similar(elevation, T)

    for i in 2:nrows-1
        for j in 2:ncols-1
            dzdx = (elevation[i, j + 1] - elevation[i, j - 1]) / (2 * cellsize)
            dzdy = (elevation[i + 1, j] - elevation[i - 1, j]) / (2 * cellsize)

            slope[i, j] = atan(sqrt(dzdx^2 + dzdy^2)) * (180 / π)
        end
    end

    return slope
end


function calculate_D(precip:Float64,elevation:::Matrix{T}, cellsize::T)
    slope = calculate_slope(elevation,cellsize)
    # function format
    return 
end

end