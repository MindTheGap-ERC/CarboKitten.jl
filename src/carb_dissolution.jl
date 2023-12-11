# based on Kaufman 2002, Geomorphology
#module carb_dissolution

export dissolution

using CarboKitten.Burgess2013
#include("climateconfig.jl") # using or import is not working, so I use include keyword

mutable struct climate
    P::Float64
    T::Float64
    pco2::Float64
end

testclimate = climate(288,1000,10^(-1.5))
const A = -0.4883+8.074*0.0001*(testclimate.T-273)
const B = -0.3241+1.6*0.0001*(testclimate.T-273)
const IA = 0.1 # ion activity
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
=#
 chemparam = Dict(
 :K1 => 10^(-356.3094-0.06091964*testclimate.T+21834.37/testclimate.T+126.8339*log10(testclimate.T)-1684915/(testclimate.T^2)),
 :K2 => 10^(-107.881-0.03252849*testclimate.T+5151.79/testclimate.T+38.92561*log10(testclimate.T)-563713.9/(testclimate.T^2)),
 :KC => 10^(-171.9065-0.077993*testclimate.T+2839.319/testclimate.T+71.595*log10(testclimate.T)),
 :KH => 10^(108.3865+0.01985076*testclimate.T-6919.53/testclimate.T-40.4515*log10(testclimate.T)+669365/(testclimate.T^2)),
 :gama_Ca => 10^(-4A*sqrt(IA)/(1+10^(-8)*B*sqrt(IA))),
 :gama_alk => 10^(-A*sqrt(IA)/(1+5.4*10^(-8)*B*sqrt(IA))),
 :alpha => 2e-6 #changeable
)

#calculate ceq and Deq, Kaufman 2002
function calculate_ceq(chemparam, cli::climate, facies::Facies)
    ceq = (cli.pco2 .* (chemparam[:K1] .* chemparam[:KC] .* chemparam[:KH]) ./(4 * chemparam[:K2] .* chemparam[:gama_Ca] .* (chemparam[:gama_alk]).^2)).^(1/3)
    Deq = cli.P .* facies.inf * 40 * 1000 * ceq ./ facies.density
    return ceq, Deq
end

# check whether the system reaches 
function dissolution(cli::climate, water_depth::Float64, chemparam, facies::Facies)
    z0 = water_depth
    I = cli.P .* facies.inf #assume vertical infiltration
    lambda = cli.P .* facies.inf ./ (chemparam[:alpha] .* facies.L)
    ceq, Deq = calculate_ceq(chemparam, cli, facies) # pass ceq Deq from the last function
    dis = (1 - exp(-z0./lambda)) > 0.8 ? Deq :  (I .* ceq ./facies.density) .* (1 - (lambda./z0)).* (1 - exp(-z0./lambda))
    return dis
end

#dissolution(testclimate,10.0,chemparam,CaProd.Facies((4, 10), (6, 10), 500, 0.8, 300, 1000, 2.73, 0.5))