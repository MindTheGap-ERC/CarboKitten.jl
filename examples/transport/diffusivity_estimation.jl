# Estimate location and width of the production peak
# using moments of the distribution, based on:
# https://github.com/MindTheGap-ERC/CarboKitten.jl/issues/110#issuecomment-2754962558
# https://gist.github.com/jhidding/1c1f0248de722fb26b356e00aa20ff5a

module DiffusivityEstimation

using CarboKitten.Export: read_slice, Header, DataSlice
using CarboKitten: MemoryOutput
using CarboKitten.Utility: in_units_of
using Unitful

# functions from Johan's gist

trapezoid(t, y) = (t[2:end] .- t[1:end-1]) .* (y[2:end] .+ y[1:end-1]) ./ 2 |> sum

moment(g, t, y) = trapezoid(t, y .* g.(t))

"""
Compute the location (μ) and width (σ) of a sediment peak from
the spatial profile `h(x)` using moments.

- `x`: position vector (with units)
- `h`: sediment thickness vector (with units)

Returns `(μ, σ)` — the center of mass and standard deviation.
"""
function peak_moments(x, h)
    s = moment(t -> 1, x, h)           # zeroth moment (total mass)
    μ = moment(t -> t, x, h ./ s)      # first moment (mean position)
    v = moment(t -> (t - μ)^2, x, h ./ s)  # second central moment (variance)
    σ = √(v)
    return (μ=μ, σ=σ)
end

"""
    peak_evolution(header::Header, data::DataSlice; use_active_layer=true)

Compute peak location and width at each saved time step
from `Header` and `DataSlice`.

Returns `(t, μ, σ)` vectors.
"""
function peak_evolution(header::Header, data::DataSlice; use_active_layer=true)
    x = header.axes.x
    t = header.axes.t[1:data.write_interval:end]

    n_t = length(t)
    Length = typeof(1.0u"m")
    μ = Vector{Length}(undef, n_t)
    σ = Vector{Length}(undef, n_t)

    if use_active_layer && data.active_layer !== nothing
        # active_layer has shape [n_facies, n_x, n_t] — sum over facies
        al = dropdims(sum(data.active_layer; dims=1); dims=1)
    else
        al = nothing
    end

    for i in 1:n_t
        h = al !== nothing ? al[:, i] : data.sediment_thickness[:, i]
        if all(h .<= 0.0u"m")
            μ[i] = NaN * u"m"
            σ[i] = NaN * u"m"
        else
            m = peak_moments(x, h)
            μ[i] = m.μ
            σ[i] = m.σ
        end
    end

    return (t=t, μ=μ, σ=σ)
end

"""
    peak_evolution(filename, group=:profile; use_active_layer=true)

Read an HDF5 file and compute peak location and width at each saved time step.

Returns `(t, μ, σ)` vectors.
"""
function peak_evolution(filename::AbstractString, group=:profile; use_active_layer=true)
    header, data = read_slice(filename, group)
    peak_evolution(header, data; use_active_layer)
end

"""
    peak_evolution(output::MemoryOutput, group=:profile; use_active_layer=true)

Compute peak location and width at each saved time step
from MemoryOutput.

Returns `(t, μ, σ)` vectors.
"""
function peak_evolution(output::MemoryOutput, group::Symbol=:profile; use_active_layer=true)
    data = output.data_slices[group]
    peak_evolution(output.header, data; use_active_layer)
end

"""
    linear_fit(x, y)

Ordinary least-squares fit of y = a + b·x.
Returns `(a, b, R²)`.
"""
function linear_fit(x, y)
    n = length(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(x .* x)
    sxy = sum(x .* y)
    syy = sum(y .* y)
    denom = n * sxx - sx^2
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    ss_res = syy - a * sy - b * sxy
    ss_tot = syy - sy^2 / n
    R² = 1 - ss_res / ss_tot
    return (a=a, b=b, R²=R²)
end

"""
    _fit_diffusivity(t, μ, σ; R²_threshold=0.99)

Internal: fit the diffusion coefficient from peak evolution data.
"""
function _fit_diffusivity(t, μ, σ; R²_threshold=0.99)

    # Filter out NaN entries (e.g. t=0 with no sediment)
    valid = .!isnan.(σ)
    t_v = t[valid]
    σ_v = σ[valid]
    μ_v = μ[valid]

    # σ² as a function of time
    σ² = σ_v .^ 2

    # Fit σ²(t) = σ₀² + 2D·t, trimming saturated tail
    n = length(t_v)
    n_fit = n
    fit = linear_fit(t_v, σ²)

    while fit.R² < R²_threshold && n_fit > 3
        n_fit -= 1
        fit = linear_fit(t_v[1:n_fit], σ²[1:n_fit])
    end

    D = fit.b / 2

    @info "Estimated diffusivity: $D"
    @info "Linear fit R² = $(fit.R²), using $(n_fit)/$(n) valid time steps"
    @info "Fit range: t = $(t_v[1]) to $(t_v[n_fit])"
    @info "Peak center at t_first: $(μ_v[1]), t_last: $(μ_v[n_fit])"
    @info "Peak width σ at t_first: $(σ_v[1]), t_last: $(σ_v[n_fit])"

    return (D=D, R²=fit.R², n_fit=n_fit, t=t, μ=μ, σ=σ)
end

"""
    estimate_diffusivity(filename, group=:profile; use_active_layer=true, R²_threshold=0.99)

Estimate the effective diffusion coefficient from an HDF5 file.
For pure diffusion, σ²(t) = σ₀² + 2Dt,
so D = slope / 2 from a linear fit of σ²(t).

The fit automatically trims late time steps where σ² saturates
(e.g. due to boundary effects) by progressively removing points
from the end until R² exceeds `R²_threshold`.

Returns a named tuple with `D`, `R²`, `n_fit`, `t`, `μ`, `σ`.
"""
function estimate_diffusivity(filename::AbstractString, group=:profile; use_active_layer=true, R²_threshold=0.99)
    t, μ, σ = peak_evolution(filename, group; use_active_layer)
    _fit_diffusivity(t, μ, σ; R²_threshold)
end

"""
    estimate_diffusivity(output::MemoryOutput, group=:profile; use_active_layer=true, R²_threshold=0.99)

Estimate the effective diffusion coefficient from MemoryOutput.
See the `AbstractString` method for details.
"""
function estimate_diffusivity(output::MemoryOutput, group::Symbol=:profile; use_active_layer=true, R²_threshold=0.99)
    t, μ, σ = peak_evolution(output, group; use_active_layer)
    _fit_diffusivity(t, μ, σ; R²_threshold)
end

end # module
