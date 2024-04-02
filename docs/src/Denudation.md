# Denudation

The denudation could be achieved by three ways: modelling physical erosion, modelling  chemical dissolutionand estimating total denudation rates based on Chlorine (Cl) isotope data.

## Physical Erosion and sediments redistribution
This method, not only considers the amount of materials that have been removed, but also how the eroded materials being distributed to the neighboring regions depending on slopes on each direction.

### Physical Erosion 
The equations used to estimate how much material could one cell provide to the lower cells is described underneath. The equation is found in [Tucker et al., 1998](https://link.springer.com/chapter/10.1007/978-1-4615-0575-4_12). We choose this equation is mainly because this equation specifically is dealing with bedorck substrates instead of loose sediments. In the equation, kv is erodibility, and the default value is 0.23 according to the paper. (1 - I_f) indicates run-off generated in one cell and slope is the slope calculated based on [ArcGis: how slope works](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm). 

$$D_{phys} = -k_v * (1 - I_f)^{1/3} |\nabla h|^{2/3}$$

``` {.julia #physical-erosion}
function physical_erosion(slope::Float64, Facies::facies)
    local kv = 0.23 #very arguable paramster
    #stencil(Float64,Reflected{2},(3,3),function(w)
    -kv .* (1-Facies.inf).^(1/3) .* slope.^(2/3)
end
```

### Redistribution of sediments

The redistribution of sediments after physical erosion is based on [van der Wiel et al., 2007](https://www.sciencedirect.com/science/article/pii/S0169555X07001341): the eroded sediments that calculated from the above equation are distributed to the neighboring 8 cells according to the slopes (defined as elevation differences/horizontal differences) towards each direction. The amount of sediments of one cell received is calculated by three functions below: 

- find the kernel to calculate redistibution co-efficient for the neighboring 8 cells depending on slopes

``` {.julia file=src/Erosion.jl}
<<physical-erosion>>
<<erosion-transport>>
```

``` {.julia #erosion-transport}
function redistribution_kernel(w::Matrix{Float64},cellsize::Float64)
    s = zeros(Float64,(3,3))
	s[1,1] = -(w[1,1] - w[2,2]) / cellsize
    s[1,2] = -(w[1,2] - w[2,2]) / cellsize / sqrt(2)
    s[1,3] = -(w[1,3] - w[2,2]) / cellsize
    s[2,1] = -(w[2,1] - w[2,2]) / cellsize / sqrt(2) 
    s[2,2] = -(w[2,2] - w[2,2]) / cellsize
    s[2,3] = -(w[2,3] - w[2,2]) / cellsize / sqrt(2)
    s[3,1] = -(w[3,1] - w[2,2]) / cellsize 
    s[3,2] = -(w[3,2] - w[2,2]) / cellsize / sqrt(2)
    s[3,3] = -(w[3,3] - w[2,2]) / cellsize

	for i in CartesianIndices(s)
		if s[i] > 0
		   continue
		else
		   s[i] = 0.0
		end
	end
	sumslope = sum(s)

	if sumslope == 0.0
	zeros(Float64,(3,3))
	else
	s./sumslope
	end
end
```

- find how much sediments would distributed to the neighboring 8 cells

``` {.julia #erosion-transport}
function mass_erosion(::Type{T},::Type{BT},slope::Matrix{Float64},n::NTuple{dim,Int}) where {T, dim, BT <: Boundary{dim}}
	m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)
	redis = zeros(Float64,(3,3,size(slope)...))
	local inf = 0.5
	for i in CartesianIndices(slope)
	     #println(i)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
			#println(Δi)
            stencil[k] = offset_value(BT, w, i, Δi)
			#println(k)
			redis[:,:,i] .= -1 .* redistribution_kernel(stencil,csz) .* physical_erosion(slope[i],inf)
        end
    end
	return redis		

end
```

- how much sediments would one cell receive in total

``` {.julia #erosion-transport}
function total_mass_redistribution(redis::Array{Float64},slope::Matrix{Float64})
	result = zeros(Float64,size(slope))
	for idx in CartesianIndices(redis)
		for i in CartesianIndices(slope)
			#println(idx)
			#println(i)
			#println(idx[1],idx[2],idx[3],idx[4],i[1],i[2])
		if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
		result[i] += redis[idx]
		end
		end
	end
	return result
end
end
```

## Chemical dissolution
The details could be found in paper by Kaufman 2002, Terra Nova. 
[link is here](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-3121.2001.00345.x)

Limestone is made of CaCO3, easily dissolved. This depends mainly on precipitation (rainfall) and temperature. The paper used equation 1 to quantify this process.

$(dh/dt) = 0.001 * kc * qi/Ai$ 

Herein dh/dt is the chemical weathering rate and the unit is in m/s.
Other parameters are defined as: qi is the discharge of water at a certain cell. Ai is the surface area od the cell. If we assume there would be no surface water on land, the qi reduces to precipitation – evaporation. Let’s set it to 400mm/y for now. Therefore equation 1 could be reduced to Equation 2.

$(dh/dt) = 0.001 * kc * I$

The parameter kc is a dimensionless parameter and should be described by equation 3:

$kc = 40 * 1000 * [Ca2+]eq/ρ$

Parameter ρ is density of calcite, and we choose 2700 kg/m3 here. The [Ca2+]eq is defined in equation 4:

$[Ca2+]eq = (PCO2 * (K1 * KC * KH) / (4 * K2 * γCa * (γHCO3)^2))^{(1/3)}$

These K1, K2, KC, KH depends on temperature. PCO2 could be between 10^(-1.5) to 10^(-3.5). 

Other parameters could be found in the following figure.

![image](https://github.com/MindTheGap-ERC/CarboKitten/assets/64159957/0f8af162-c508-4873-b5ab-25d39a955f87)

However, the above discussion is true only if the percolated fluid is saturated (in terms of Ca) when leaving the platform. In some cases, when the fluid is not saturated, the dissolvbed amount is lower than the scenario described above.

The following Paper describes this:
[click here](https://www.sciencedirect.com/science/article/pii/S0169555X08004133#bib22)
and/or
[click here](https://www.sciencedirect.com/science/article/pii/S001670370602240X)

Ideally, a reactive transport model should be most accurate but it needs more calculation resources. So herein, the author just suggested the dissolution rates of rocks depend on the depth. This makes sense as the deeper you go the more concentrated you are. Also, this does not consider diffusion in this chapter.

The dissolution rate of carbonate follows linear rate laws of:
$F = α (ceq-c(z))$

The rate law is a common expression way to describe the kinetics of certain chemical reactions and [click hier for more info](https://chem.libretexts.org/Bookshelves/General_Chemistry/Chemistry_1e_(OpenSTAX)/12%3A_Kinetics/12.3%3A_Rate_Laws)

F is the dissolution rate, α is constant (kinetic co-efficient), ceq is the concentration in fluid when equilibrium is reached (i.e., no more dissolution, which is [Ca2+]eq in Chapter 1), c(z) is the current concentrationion at depth of z in the fluid. This equation then expands to 

$I *dc = α (ceq-c(z)) * L dz$

This equation indicates that the concentration increases in the infiltrated water equals the dissolution of rocks in the thickness of 'dz'. 'L' is the specific length of fractures/porosities (units: m/m^2, we can try 100 at the first place). I.e., this term defines the relative reactive surface of the subsurface rocks, or how much surface is actually dissolving? This term is difficult to determine. 'I' is infiltration, but slightly different as chapter 1: this I is the I in each rain event according to the paper. We certainly do not gonna know how this parameter works, so we just set it the same as in chapter 1?

However, to solve this equation we still need to know c(z).

If assuming the initial percolating water has c(0) = 0, then we could get the following equation (as the c is related to depth):

$c(z) = ceq * (1 - e^{(-z/λ)})$

Herein, λ = I/αL. 

Therefore the $Daverage = (I * ceq/ρ) * (1 – (λ/z0) * (1 – e^{(-z0/λ)}))$

α used in this article is α = 2·10^{−6} or 3.5·10^{−7} cm/s (for temp at 298K). This is indeed a controversial parameter TBH. We can try different values and see what happens.

## Emperical denudation
Cl isotopes are an emerging tool to decipher the denudation rates (chemical dissolution + physical erosion) in carbonate-dominated area.

Research based on the karst region and carbonate platform terrace suggested that the denudation rates are mainly controlled by precipitation and slopes, although the debates about which factor is more important is still ongoing ([Ryb et al., 2014](https://www.sciencedirect.com/science/article/pii/S1871101420300248#sec4), [Thomas et al., 2018](https://www.sciencedirect.com/science/article/pii/S0169555X18301740)). In general, the precipitation mainly controls the chemical dissolution while the slopes mainly controls the physical ersions. In addition, the type of carbonates may also play an important role ([Kirklec et al., 2022](https://www.sciencedirect.com/science/article/pii/S0169555X22002513)), but given this feature is studied poorly so we will ditch it for now. We have checked and compiled the denudation rates (mm/kyr) along with precipitation and slopes serve as a starting point to create a function relates denudation rates (mm/kyr) to precipitation and slopes. The compiled data could be found in OSFdatabase. This is an empirical relationship and have a relatively large uncertainty in terms of fitting. 
<img width="433" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/29da37f3-0e14-479f-8ed3-a8b5adc052a5">
*Fig 1. The relationship between MAP (mean precipitation per year, mm/y) and denudation rates (mm/ky)*
<img width="370" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/1bd9e2ae-d4d7-4611-a8b8-83459222022f">
*Fig 2. The relationship between the curvature (i.e., slope or height) and the denudation rates (mm/ky)*

We can see that both the slope and precipitation could increase the denudation rates, and reaches a 'steady state' after a certain point.

In terms of impletementation, I used the function form of D = P * S, where D means denudation rates, P means effects of precipitation while S means effects of Slope. By doing so, we can consider both effects. Such formula structure is similar to RUSLE model, a widely used LEM (e.g., (Thapa et al., 2020)[https://environmentalsystemsresearch.springeropen.com/articles/10.1186/s40068-020-00177-2]).

The codes are attached:

## How to use?
In CarboKitten, you could choose which type of the three you would like to attempt. To do this you could simply change the 'erosino_type' in the input. 

Example: in file examples, you could find caps_miller. This file use [Miller's curve](https://www.science.org/doi/full/10.1126/science.1116412) as sealevel curve input. You could simply try different erosion types by changing the 'erosion type':
- 0 means no erosion
- 1 means chemical dissolution
- 2 means physical erosion and sediments redistribution
- 3 means total denudation calculated based on emperical relationship by Cl isotope observations.
