# Denudation

Denudation could be achieved by three ways: modelling physical erosion, modelling chemical dissolution and estimating total denudation rates based on chlorine (Cl) isotope data.

## Physical erosion and sediment redistribution
This method not only considers the amount of materials that have been removed, but also how the eroded materials being distributed to the neighboring regions depending on slopes on each direction.

### Physical erosion 
The equations used to estimate how much material could one cell provide to the lower cells is described underneath. The equation is found in [tucker_channel-hillslope_2001](@cite). We choose this equation mainly because it specifically deals with bedrock substrates instead of loose sediments. In the equation, $k_v$ is erodibility, and the default value is 0.23 according to the paper. $(1 - I_f)$ indicates run-off generated in one cell and slope is the slope calculated based on [ArcGis: how slope works](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm). Note that the algorithms to calculate slope does not work on depressions. 

$$D_{phys} = -k_v * (1 - I_f)^{1/3} |\nabla h|^{2/3}$$

``` {.julia #physical-erosion}
function physical_erosion(slope::Float64, Facies::facies)
    local kv = 0.23 #very arguable paramster
    #stencil(Float64,Reflected{2},(3,3),function(w)
    -kv .* (1-Facies.inf).^(1/3) .* slope.^(2/3)
end
```

### Redistribution of sediments

The redistribution of sediments after physical erosion is based on [van_de_wiel_embedding_2007](@cite): the eroded sediments calculated from the above equation are distributed to the neighboring 8 cells according to the slopes (defined as elevation differences/horizontal differences) towards each direction. The amount of sediments of one cell received is calculated by three functions below: 

#### Find the kernel to calculate redistibution co-efficient for the neighboring 8 cells depending on slopes

``` {.julia}
module Erosion

<<physical-erosion>>
<<erosion-transport>>
end  # module
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

#### Find out how much sediments would distributed to the neighboring 8 cells

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

#### How much sediment would one cell receive in total

``` {.julia #erosion-transport}
function total_mass_redistribution(redis::Array{Float64},slope::Any,::Type{BT}) where {BT <: Boundary}
	mass = zeros(Float64,size(slope))
    for i in CartesianIndices(slope)
        for idx in CartesianIndices(redis)
            if offset_index(BT, size(slope), CartesianIndex(idx[3],idx[4]), CartesianIndex(idx[1]-2,idx[2]-2)) == i
            @show i
			mass[i] += redis[idx]
            end
		end
	end
	return mass
end
```

## Chemical dissolution

The details could be found in paper by [Kaufmann2001](@cite).

Limestone is made of $CaCO_3$, easily dissolved. This depends mainly on precipitation (rainfall) and temperature. The paper used equation 1 to quantify this process.

$(dh/dt) = 0.001 * κ_c * q_i/A_i$ 

Herein $dh/dt$ is the chemical weathering rate and the unit is in m/s.
Other parameters are defined as: $q_i$ is the discharge of water at a certain cell. $A_i$ is the surface area of the cell. If we assume there would be no surface water on land, $q_i$ reduces to precipitation – evaporation. Let’s set it to 400 mm/y for now. Therefore equation 1 could be reduced to Equation 2.

$(dh/dt) = 0.001 * κ_c * I$

Where $I$ is runoff (mm/y?). The parameter $κ_c$ is dimensionless and should be described by equation 3:

$kc = 40 * 1000 * [Ca^{2+}]_{eq}/ρ$

Parameter ρ is the density of calcite, and we choose 2700 $kg/m^3$ here. $[Ca^{2+}]_{eq}$ is defined in equation 4:

$[Ca^{2+}]_{eq} = (PCO_2 * (K_1 * K_C * K_H) / (4 * K_2 * γCa * (γHCO_3)^2))^{(1/3)}$

Mass balance coefficients $K_1$, $K_2$, $K_C$, $K_H$ depend on temperature. $PCO_2$ is assumed to be between $10^{-1.5}$ to $10^{-3.5}$ (units?). 

Other parameters could be found in the following table by [Kaufmann2001](@cite).

| Parameter                 | Description              | Unit                  | Value                                                                      |
|-----------------|-----------------|-----------------|----------------------|
| *T*                       | Absolute temperature     | \[°K\]                | Tc + 273.16                                                                |
| *I*                       | Ion activity             | \[-\]                 | 0.1                                                                        |
| $A^*$                     | Debye-Hückel coefficient | \[-\]                 | $-0.4883 + 8.074 \times 10^{-4}T_c$                                             |
| $B^*$                     | Debye-Hückel coefficient | \[-\]                 | $-0.3241 + 1.600 \times 10^{-4}T_c$                                             |
| $log \gamma Ca\dagger$    | Activity coefficient     | \[-\]                 | $-4A\sqrt{I}(1 + 5.0 x 10^{-8}B\sqrt{I})$                                  |
| $log \gamma HCO_3\dagger$ | Activity coefficient     | \[-\]                 | $-1A\sqrt{I}/(1 + 5.4 x 10^{-8}B\sqrt{I})$                                 |
| $log K_1\ddagger$         | Mass balance coefficient | \[ $mol L^{-1}$ \]         | $-356.3094 - 0.06091964T + 21834.37/T + 126.8339logT - 1684915/T^2$ |
| $log K_2\ddagger$         | Mass balance coefficient | \[ $mol L^{-1}$ \]         | $-107.8871 - 0.03252849T + 5151.79/T + 38.92561logT - 563713.9/T^2$ |
| $log K_c\ddagger$         | Mass balance coefficient | \[ $mol^2 L^{-2}$ \]      | $-171.9065 - 0.077993T + 2839.319/T + 71.595logT$                      |
| $log K_H\ddagger$         | Mass balance coefficient | \[ $mol L^{-1} atm^{-1}$ \] | $108.3865 + 0.01985076T - 6919.53/T - 40.45154logT + 669365/T^2$    |

However, the above discussion is true only if the percolated fluid is saturated (in terms of Ca) when leaving the platform. In some cases, when the fluid is not saturated, the dissolved amount is lower than the scenario described above.

The following articles describe this: [gabrovsek_concepts_2009](@cite) and [kaufmann_calcite_2007](@cite)

Ideally, a reactive transport model should be accurate, but that needs more computation resources. So herein, the author just suggested the dissolution rates of rocks depend on the depth. This makes sense, as the deeper the solution penetrates, the more concentrated it becomes. Also, this does not consider diffusion in this chapter.

The dissolution rate of carbonate follows linear rate laws of:
$F = α (c_{eq}-c(z))$

The rate law is a common expression way to describe the kinetics of certain chemical reactions [see Rate Laws](https://chem.libretexts.org/Bookshelves/General_Chemistry/Chemistry_1e_(OpenSTAX)/12%3A_Kinetics/12.3%3A_Rate_Laws)

$F$ is the dissolution rate, α is constant (kinetic co-efficient), $c_{eq}$ is the concentration in fluid when equilibrium is reached (i.e., no more dissolution, which is $[Ca^{2+}]_{eq}$ in Chapter 1), $c(z)$ is the current concentrationion at depth $z$ in the fluid. This equation then expands to 

$I *dc = α (c_{eq}-c(z)) * L dz$

This equation indicates that the concentration increase in the infiltrated water equals the dissolution of rocks in the thickness of $dz$. $L$ is the specific length of fractures/porosities (units: $m/m^2$, we can try 100 at the first place). I.e., this term defines the relative reactive surface of the subsurface rocks, or how much surface is actually dissolving. This term is difficult to determine. $I$ is infiltration, but slightly different as chapter 1: this $I$ is the $I$ in each rain event according to the paper. We certainly do not gonna know how this parameter works, so we just set it the same as in chapter 1?

However, to solve this equation we still need to know $c(z)$.

If assuming the initial percolating water has $c(0) = 0$, then we could get the following equation (as $c$ is related to depth):

$c(z) = c_{eq} * (1 - e^{(-z/λ)})$

Herein, λ = I/αL. 

Therefore $D_{average} = (I * c_{eq}/ρ) * (1 – (λ/z0) * (1 – e^{(-z0/λ)}))$

α used in this article is $α = 2·10^{−6}$ or $3.5·10^{−7}$ cm/s (for temp at 298K). This is indeed a controversial parameter TBH. We can try different values and see what happens.

## Emperical denudation
Cl isotopes are an emerging tool to decipher the denudation rates (chemical dissolution + physical erosion) in carbonate-dominated area.

Research based on the karst region and carbonate platform terrace suggested that the denudation rates are mainly controlled by precipitation and slopes, although the debates about which factor is more important is still ongoing ([yang_combined_2020](@cite), [thomas_limited_2018](@cite)). In general, the precipitation mainly controls the chemical dissolution while the slopes mainly controls the physical ersions. In addition, the type of carbonates may also play an important role ([krklec_long-term_2022](@cite)), but given this feature is studied poorly so we will ditch it for now. We have checked and compiled the denudation rates (mm/kyr) along with precipitation and slopes serve as a starting point to create a function relates denudation rates (mm/kyr) to precipitation and slopes. The compiled data could be found in OSFdatabase. This is an empirical relationship and have a relatively large uncertainty in terms of fitting.

<img width="433" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/29da37f3-0e14-479f-8ed3-a8b5adc052a5">

*Fig 1. The relationship between MAP (mean precipitation per year, mm/y) and denudation rates (mm/ky)*

<img width="370" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/1bd9e2ae-d4d7-4611-a8b8-83459222022f">

*Fig 2. The relationship between the curvature (i.e., slope or height) and the denudation rates (mm/ky)*

We can see that both the slope and precipitation could increase the denudation rates, and reaches a 'steady state' after a certain point.

In terms of implementation, I used the function form of $D = P * S$, where $D$ means denudation rates, $P$ means effects of precipitation while $S$ means effects of Slope. By doing so, we can consider both effects. Such formula structure is similar to RUSLE model, a widely used LEM (e.g., (thapa_spatial_2020)[@cite]).

## How to use?
In CarboKitten, you could choose which type of the three you would like to attempt. To do this you could simply change the `erosion_type` in the input. 

Example: in examples, you find `caps_miller_diss.jl`, `caps_miller_emp.jl`, `caps_miller_phys.jl`, for chemical dissolution, empirical denudation or physical denudation, respectively. This file uses the [Mmiller_phanerozoic_2005] (@cite) cure as sea level curve input. You could try different erosion types by changing the `erosion_type`:
- `NoDenudation` means no erosion
- `Dissolution` means chemical dissolution
- `PhysicalErosionParam` means physical erosion and sediments redistribution
- `EmpericalDenudationParam` means total denudation calculated based on emperical relationship by Cl isotope observations.
