module AriealE

using ..Production: w
using ..Transport: elevation

struct para
    silt::Float16
    sand::Float16
    clay::Float16
    precip::Float64
end

function slopecalculation (elevation)
    slope = -elevation/100
end
# silt, sand and clay are proportions

function physicalE
    A = R * K * L * S 
    R = 1.735 * 10 ^(1.5 * log(para.precip)-0.8188)
    K = (0.2 + 0.3exp(-0.0256*para.sand*(1-para.silt))) * (para.silt/(para.silt+para.clay))^0.3 * (1-(0.7*(1-para.sand)/((1-para.sand)+exp(-5.51+22.9*(1-para,sand)))))
    S = 10.8 * sin(slope) + 0.03
    beta = sin(slope)/(3*sin(slope)^0.8 + 0.56)
    alpha = (beta/(beta+1))
    L = (100/22.3)^alpha
end

end

