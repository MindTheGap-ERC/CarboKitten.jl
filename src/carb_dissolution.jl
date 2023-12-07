# based on Kaufman 2002, Geomorphology
#module carb_dissolution

export dissolution
using CarboKitten.Burgess2013
import climateconfig: climate, model1

const A = -0.4883+8.074*0.0001*(climate.T-273)
const B = -0.3241+1.6*0.0001*(climate.T-273)
const IA = 0.1 # ion activity

# chemical basic parameters
struct chemparam
    K1::Float64
    K2::Float64
    KC::Float64
    KH::Float64
    gama_Ca::Float64
    gama_alk::Float64
    a::Float64
end

const modeldis = [chemparam(10^(-356.3094-0.06091964*climate.T+21834.37/climate.T+126.8339*log(climate.T)-1684915/(climate.T^2)), 
10^(-107.881-0.03252849*climate.T+5151.79/climate.T+38.92561*log(climate.T)-5637.13/(climate.T^2)), 
10^(-171.9065-0.077993*climate.T+2839.319*climate1.T+71.595*log(climate.T)), 
10^(108.3865-0.01985076*climate.T-6919.53/climate.T+40.4515*log(climate.T)-669365/(climate.T^2)),
(-4A*sqrt(IA)*(1+5e-8*B*sqrt(IA))),
(-A*sqrt(IA)*(1+5.4e-8*B*sqrt(IA))))]

#calculate ceq and Deq, Kaufman 2002
function calculate_ceq(param::chemparam, cli::climate, facies::Facies)
    ceq = (cli.pco2 .* (param.K1 .* param.KC .* param.KH) ./(4 * param.K2 .* param.gama_Ca .* (param.gama_alk).^2)).^(1/3)
    Deq = 0.001 * cli.P .* facies.inf * 40 * 1000 * ceq ./ facies.density
    return ceq, Deq
end

# check whether the system reaches 
function dissolution(cli::Climate, water_depth::Float64, param::chemparam,  facies::Facies)
    z0 = water_depth
    I = cli.P .* facies.inf #assume vertical infiltration
    lambda = cli.P .* facies.inf ./ (param.a .* facies.L)
    ceq, Deq = calculate_ceq(param, cli, facies) # pass ceq Deq from the last function
    dis = (1 - exp(-z0./lambda)) > 0.8 ? Deq :  (I .* ceq ./facies.density) .* (1 - (lambda./z0)).* (1 - exp(-z0./lambda))
    return dis
end

#end