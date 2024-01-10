# based on Kaufman 2002, Geomorphology
module CarbDissolution

export dissolution

using ..Burgess2013

#include("climateconfig.jl") # using or import is not working, so I use include keyword
#=
mutable struct climate
    P::Float64
    T::Float64
    pco2::Float64
end
=#
#testclimate = climate(1000,288,10^(-1.5))

# Kaufmann 2002, Table 2
function karst_denudation_parameters(temp::Float64)
    A = -0.4883+8.074*0.0001*(temp-273)
    B = -0.3241+1.6*0.0001*(temp-273)
    IA = 0.1 # ion activity

    ( K1 = 10^(-356.3094-0.06091964*temp+21834.37/temp+126.8339*log10(temp)-1684915/(temp^2))
    , K2 = 10^(-107.881-0.03252849*temp+5151.79/temp+38.92561*log10(temp)-563713.9/(temp^2))
    , KC = 10^(-171.9065-0.077993*temp+2839.319/temp+71.595*log10(temp))
    , KH = 10^(108.3865+0.01985076*temp-6919.53/temp-40.4515*log10(temp)+669365/(temp^2))
    , activity_Ca = 10^(-4A*sqrt(IA)/(1+10^(-8)*B*sqrt(IA)))
    , activity_Alk = 10^(-A*sqrt(IA)/(1+5.4*10^(-8)*B*sqrt(IA)))
    )
end

#=
const T = 288
const P = 1000
const pco2 = 10^(-1.5)
const A = -0.4883+8.074*0.0001*(T-273)
const B = -0.3241+1.6*0.0001*(T-273)
const IA = 0.1 # ion activity
const a = 
=#
#= chemical basic parameters
struct chemparam
    K1::Float64
    K2::Float64
    KC::Float64
    KH::Float64
    gama_Ca::Float64
    gama_alk::Float64
    a::Float64
end

 chemparam = Dict(
 :K1 => 10^(-356.3094-0.06091964*testclimate.T+21834.37/testclimate.T+126.8339*log10(testclimate.T)-1684915/(testclimate.T^2)),
 :K2 => 10^(-107.881-0.03252849*testclimate.T+5151.79/testclimate.T+38.92561*log10(testclimate.T)-563713.9/(testclimate.T^2)),
 :KC => 10^(-171.9065-0.077993*testclimate.T+2839.319/testclimate.T+71.595*log10(testclimate.T)),
 :KH => 10^(108.3865+0.01985076*testclimate.T-6919.53/testclimate.T-40.4515*log10(testclimate.T)+669365/(testclimate.T^2)),
 :gama_Ca => 10^(-4A*sqrt(IA)/(1+10^(-8)*B*sqrt(IA))),
 :gama_alk => 10^(-A*sqrt(IA)/(1+5.4*10^(-8)*B*sqrt(IA))),
 :alpha => 2e-6 #changeable
)
=#
#calculate ceq and Deq, Kaufman 2002
function equilibrium(temp::Float64, pco2::Float64, precip::Float64, facies::Facies)
    p = karst_denudation_parameters(temp)
    eq_c = (pco2 .* (p.K1 * p.KC * p.KH) ./ (4 * p.K2 * p.activity_Ca * (p.activity_Alk)^2)).^(1/3)
    eq_d = precip .* facies.inf * 40 * 1000 * eq_c ./ facies.density
    (concentration = eq_c, denudation = eq_d)
end

# check whether the system reaches equilibrium
function dissolution(temp::Float64,precip::Float64, alpha::Float64, pco2::Float64,water_depth::Float64, facies::Facies)
        z0 = -water_depth
        I = precip .* facies.inf #assume vertical infiltration
        λ = precip .* facies.inf ./ (alpha .* facies.L)
        eq = equilibrium(temp,pco2,precip,facies) # pass ceq Deq from the last function
        if (1 - exp(-z0./λ)) > 0.8
            eq.denudation 
        else
            (I .* eq.concentration ./facies.density) .* (1 - (λ./z0).* (1 - exp(-z0./λ)))
        end
end

#dissolution(testclimate,10.0,chemparam,CaProd.Facies((4, 10), (6, 10), 500, 0.8, 300, 1000, 2.73, 0.5))
end