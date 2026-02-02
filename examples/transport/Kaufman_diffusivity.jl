# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/Kaufman_diffusivity.jl>>[init]
module KaufmanCommon

using Unitful
using CarboKitten
using CarboKitten.Models.ALCAP

export KAUFMAN_C0, KAUFMAN_C1, BACKGROUND_DIFFUSION
export kaufman_diffusivity, print_diffusivity_profile

# note to self: I changed the unit from m^2/yr to m/yr to match CK's units - this needs to be clarified
const KAUFMAN_C0 = 5000.0u"m/yr" # Surface diffusivity for carbonates
const KAUFMAN_C1 = 0.05u"m^-1" # Depth decay constant
const BACKGROUND_DIFFUSION = 2.0u"m^2/yr" 

kaufman_diffusivity(depth) = KAUFMAN_C0 * exp(-KAUFMAN_C1 * depth)

function print_diffusivity_profile()
    D_surface = kaufman_diffusivity(0.0u"m")

    for depth in [0, 2, 5, 10, 20, 30, 50, 100]u"m"
        D = kaufman_diffusivity(depth)
        ratio = D / D_surface |> NoUnits
        println("  $depth   │  $D  │  $ratio")
    end
end

end 
# ~/~ end
