module WaveField

using Unitful
using Unitful: ustrip, @u_str, uconvert
using GeometryBasics: Vec2

export WaveComponent, WaveModel, wave_physics_at_cell

const ρ = 1025.0u"kg/m^3"
const g = 9.81u"m/s^2"
const Myr_to_s = 3.1536e13
const γ_breaker = 0.78

struct WaveComponent
    A0          # Amplitude (m)
    T           # Period (s)
    direction   # Radians
    phase       # Radians
    mu_base     # Base attenuation (1/m)
end
struct WaveModel
    components::Vector{WaveComponent}
    efficiency
end

dir(θ) = Vec2(cos(θ), sin(θ))


function _wave_velocity_at_depth(model::WaveModel, pos, t, h, x_fetch)

    v = Vec2(0.0u"m/Myr", 0.0u"m/Myr")
    h_safe = max(h, 0.01u"m")

    for comp in model.components

        # --- wavelength from period ---
        λ = wavelength_from_period(comp.T, h_safe)

        # --- attenuation ---
        μ = comp.mu_base * (1.0u"m" / h_safe)
        A_local = comp.A0 * exp(-ustrip(μ * x_fetch))

        H_local = 2 * A_local
        H_max = γ_breaker * h_safe
        H = min(H_local, H_max)

        # --- orbital velocity ---
        T = comp.T
        denom = sinh(ustrip(2π * h_safe / λ))
        Ub = (π * H) / (T * denom)

        U_geo = Ub * Myr_to_s * model.efficiency

        ω = 2π / ustrip(u"s", T)
        phase = 2π * ustrip(pos / λ) - ustrip(ω * t) + comp.phase

        v += dir(comp.direction) .* (U_geo * sin(phase))
    end

    return v
end

const g_val = ustrip(u"m/s^2", g)  # m/s^2 (unitless for solver)


function wavelength_from_period(T, h;
                                tol=1e-10,
                                maxiter=100)

    T_s = ustrip(u"s", T)
    h_m = ustrip(u"m", h)

    ω = 2π / T_s

    # Deep water initial guess
    k = ω^2 / g_val

    for _ in 1:maxiter
        f  = g_val*k*tanh(k*h_m) - ω^2
        df = g_val*tanh(k*h_m) +
             g_val*k*h_m*(1 / cosh(k*h_m))^2
        k_new = k - f/df
        abs(k_new - k) < tol && break
        k = k_new
    end

    λ = 2π / k
    return λ * u"m"
end


function wave_physics_at_cell(model::WaveModel, pos, t, h, x_fetch; dh=0.05u"m")

    v = _wave_velocity_at_depth(model, pos, t, h, x_fetch)
    v_plus = _wave_velocity_at_depth(model, pos, t, h + dh, x_fetch)
    s = (v_plus - v) / dh

    total_p = 0.0u"kW/m"
    h_safe = max(h, 0.01u"m")

    for comp in model.components

    μ = comp.mu_base * (1.0u"m" / h_safe)
    A_local = comp.A0 * exp(-ustrip(μ * x_fetch))

    H_local = 2 * A_local
    H_max = γ_breaker * h_safe
    H = min(H_local, H_max)

    λ = wavelength_from_period(comp.T, h_safe)
    k = 2π / λ
    ω = 2π / comp.T

    c = ω / k
    n = 0.5 * (1 + (2k*h_safe) / sinh(2k*h_safe))
    c_g = c * n

    P = (0.125 * ρ * g * H^2) * c_g
    total_p += uconvert(u"kW/m", P)
end

    return (v, s, total_p)
end

end