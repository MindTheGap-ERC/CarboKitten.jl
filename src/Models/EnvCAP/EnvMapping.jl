# ~/~ begin <<docs/src/models/envcap.md#src/Models/EnvCAP/EnvMapping.jl>>[init]
module EnvMapping

using Unitful

export dominant_env, env_to_factory_prior, apply_prior_bias!

# ~/~ begin <<docs/src/models/envcap.md#dominant-env>>[init]
function dominant_env(deposition::AbstractArray{T,4}) where T
    total = dropdims(sum(deposition; dims=4), dims=4)  # (n_envs, nx, ny)
    n_envs, nx, ny = size(total)
    result = zeros(Int, nx, ny)
    for ix in 1:nx, iy in 1:ny
        best_env = 0
        best_val = zero(T)
        for e in 1:n_envs
            v = total[e, ix, iy]
            if v > best_val
                best_val = v
                best_env = e
            end
        end
        result[ix, iy] = best_env
    end
    return result
end
# ~/~ end
# ~/~ begin <<docs/src/models/envcap.md#env-to-factory-prior>>[init]
function env_to_factory_prior(
        env_field::Matrix{Int},
        mapping::Matrix{Float64})::Array{Float64,3}

    n_envs, n_factories = size(mapping)
    nx, ny = size(env_field)

    normalised = similar(mapping)
    for e in 1:n_envs
        s = sum(mapping[e, :])
        normalised[e, :] .= s > 0.0 ? mapping[e, :] ./ s : fill(1.0/n_factories, n_factories)
    end

    prior = zeros(Float64, n_factories, nx, ny)
    for ix in 1:nx, iy in 1:ny
        e = env_field[ix, iy]
        if e == 0 || e > n_envs
            prior[:, ix, iy] .= 1.0 / n_factories
        else
            prior[:, ix, iy] .= normalised[e, :]
        end
    end
    return prior
end
# ~/~ end
# ~/~ begin <<docs/src/models/envcap.md#apply-prior-bias>>[init]
function apply_prior_bias!(
        ca::Matrix{Int},
        factory_prior::Array{Float64,3},
        ca_refinement::Float64,
        rng)

    ca_refinement == 0.0 && return

    for ix in axes(ca, 1), iy in axes(ca, 2)
        f = ca[ix, iy]
        f == 0 && continue
        p_prior = factory_prior[f, ix, iy]
        p_kill  = ca_refinement * (1.0 - p_prior)
        if rand(rng) < p_kill
            ca[ix, iy] = 0
        end
    end
end
# ~/~ end

end
# ~/~ end
