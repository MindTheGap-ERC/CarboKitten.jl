# ~/~ begin <<docs/src/components/wavefield.md#src/Components/WaveField.jl>>[init]
module WaveField

using Unitful
using GeometryBasics

export AiryWaveComponent, AiryWaveField, energy_flux, group_velocity

# =============================================================================
# Types
# =============================================================================

"""
    AiryWaveComponent(; amplitude, period, direction, phase=0.0, attenuation=0.0u"m^-1")

A single monochromatic wave component:

- `amplitude` — wave height H (crest-to-trough) in meters.
- `period` — wave period T in seconds.
- `direction` — propagation angle θ in radians (0 = +x, π/2 = +y).
- `phase` — phase offset φ₀ in radians.
- `attenuation` — spatial decay coefficient (e.g. for bottom friction).
"""
@kwdef struct AiryWaveComponent
    amplitude::typeof(1.0u"m")       = 1.0u"m"
    period::typeof(1.0u"s")          = 10.0u"s"
    direction::Float64               = 0.0
    phase::Float64                   = 0.0
    attenuation::typeof(1.0u"m^-1")  = 0.0u"m^-1"
end

"""
    AiryWaveField(; components, breaker_index=0.78, gravity=9.81u"m/s^2", density=1025.0u"kg/m^3")

A superposition of `AiryWaveComponent`s. Callable: `wf(water_depth)` returns
`(Vec2{velocity}, Vec2{shear})` compatible with the existing
`Facies.wave_velocity` interface in `Components.ActiveLayer`.

Existing lambda-based wave velocities (including the default no-transport
closure and `v_const`) continue to work — this type is an alternative, not
a replacement.

- `breaker_index` — γ_b for depth-limited breaking (default 0.78).
- `gravity` — gravitational acceleration.
- `density` — water density (for energy flux diagnostics).
"""
@kwdef struct AiryWaveField
    components::Vector{AiryWaveComponent} = AiryWaveComponent[]
    breaker_index::Float64               = 0.78
    gravity::typeof(1.0u"m/s^2")         = 9.81u"m/s^2"
    density::typeof(1.0u"kg/m^3")        = 1025.0u"kg/m^3"
end

# =============================================================================
# Dispersion relation
# =============================================================================

"""
    solve_dispersion(ω, h, g; tol, maxiter) -> k

Solve ω² = g k tanh(k h) for wave number k using Newton–Raphson iteration,
starting from the deep-water approximation k₀ = ω²/g. All arguments are
plain Float64 in SI units.
"""
function solve_dispersion(ω::Float64, h::Float64, g::Float64;
                          tol::Float64=1e-8, maxiter::Int=50)
    h <= 0.0 && return 0.0
    k = ω^2 / g                        # deep-water initial guess
    for _ in 1:maxiter
        kh = k * h
        th = tanh(kh)
        f  = ω^2 - g * k * th
        df = -g * (th + k * h * (1.0 - th^2))
        abs(df) < 1e-30 && break
        dk = -f / df
        k += dk
        abs(dk) < tol * abs(k) && break
    end
    return max(k, 0.0)
end

# =============================================================================
# Per-component orbital velocity (all SI, no units)
# =============================================================================

"""
    _orbital_velocity(H, T, k, h, γ_b) -> (u_b, H_eff)

Near-bed orbital velocity amplitude and depth-limited wave height.
All arguments in SI (m, s), returns (m/s, m).
"""
function _orbital_velocity(H::Float64, T::Float64, k::Float64, h::Float64, γ_b::Float64)
    h <= 0.0 && return (0.0, 0.0)

    # Depth-limited breaking
    H_eff = min(H, γ_b * h)

    # u_b = π H / (T sinh(kh))
    kh = k * h
    sh = sinh(max(kh, 1e-10))
    u_b = π * H_eff / (T * sh)

    return (u_b, H_eff)
end

# =============================================================================
# Velocity vector at a given depth (SI, no units)
# =============================================================================

function _velocity_si(wf::AiryWaveField, h::Float64)
    g  = ustrip(u"m/s^2", wf.gravity)
    γ_b = wf.breaker_index
    vx = 0.0
    vy = 0.0

    for comp in wf.components
        H  = ustrip(u"m", comp.amplitude)
        T  = ustrip(u"s", comp.period)
        ω  = 2π / T
        k  = solve_dispersion(ω, h, g)
        u_b, _ = _orbital_velocity(H, T, k, h, γ_b)

        vx += u_b * cos(comp.direction)
        vy += u_b * sin(comp.direction)
    end

    return (vx, vy)
end

# =============================================================================
# Callable interface: (::AiryWaveField)(water_depth) -> (Vec2, Vec2)
# =============================================================================

"""
    (wf::AiryWaveField)(water_depth)

Evaluate the wave field at `water_depth`. Returns
`(velocity::Vec2, shear::Vec2)` in units of m/s and 1/s respectively,
compatible with the `Facies.wave_velocity` interface.

Velocity is the vector sum of near-bed orbital velocity amplitudes directed
along each component's propagation direction. Shear is d(velocity)/d(depth),
computed by central finite differences.
"""
function (wf::AiryWaveField)(water_depth)
    h = ustrip(u"m", water_depth)

    vx, vy = _velocity_si(wf, h)

    # Shear = dv/dh via central differences
    δ = max(abs(h) * 1e-4, 1e-4)
    vx_p, vy_p = _velocity_si(wf, h + δ)
    vx_m, vy_m = _velocity_si(wf, max(0.0, h - δ))
    sx = (vx_p - vx_m) / (2δ)
    sy = (vy_p - vy_m) / (2δ)

    return (Vec2(vx * u"m/s", vy * u"m/s"),
            Vec2(sx * u"s^-1", sy * u"s^-1"))
end

# =============================================================================
# Diagnostics: group velocity and energy flux
# =============================================================================

"""
    group_velocity(comp, h, g) -> c_g  [m/s]

Group velocity from linear theory: c_g = c · n where
n = ½(1 + 2kh / sinh(2kh)) and c = ω/k.
"""
function group_velocity(comp::AiryWaveComponent, h_m::Float64, g::Float64)
    T  = ustrip(u"s", comp.period)
    ω  = 2π / T
    k  = solve_dispersion(ω, h_m, g)
    k <= 0.0 && return 0.0
    kh = k * h_m
    c  = ω / k
    n  = 0.5 * (1.0 + 2kh / sinh(max(2kh, 1e-10)))
    return c * n
end

"""
    energy_flux(wf::AiryWaveField, water_depth) -> E  [W/m]

Total wave energy flux per unit crest length from all components:
E = Σ (1/8) ρ g H_eff² c_g
"""
function energy_flux(wf::AiryWaveField, water_depth)
    h  = ustrip(u"m", water_depth)
    g  = ustrip(u"m/s^2", wf.gravity)
    ρ  = ustrip(u"kg/m^3", wf.density)
    γ_b = wf.breaker_index

    E = 0.0
    for comp in wf.components
        H  = ustrip(u"m", comp.amplitude)
        T  = ustrip(u"s", comp.period)
        ω  = 2π / T
        k  = solve_dispersion(ω, h, g)
        _, H_eff = _orbital_velocity(H, T, k, h, γ_b)
        c_g = group_velocity(comp, h, g)
        E += (1.0 / 8.0) * ρ * g * H_eff^2 * c_g
    end
    return E * u"W/m"
end

end
# ~/~ end
