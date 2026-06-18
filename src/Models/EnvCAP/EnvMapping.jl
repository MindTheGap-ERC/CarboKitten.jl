# ~/~ begin <<docs/src/models/envcap.md#src/Models/EnvCAP/EnvMapping.jl>>[init]
module EnvMapping

using Unitful
using Unitful: NoUnits

export dominant_env_block, env_to_factory_prior_block, apply_prior_bias!

# ~/~ begin <<docs/src/models/envcap.md#dominant-env>>[init]
function dominant_env_block(deposition::AbstractArray{T,4}; dz = 0.5, min_nz = 1) where T
    n_envs, nx, ny, nt = size(deposition)

    numeric_dep = try
        ustrip.(u"m", deposition)
    catch
        Float64.(deposition)
    end

    total_thickness = dropdims(sum(numeric_dep; dims = (1, 4)), dims = (1, 4))
    nz = max(min_nz, 1, ceil(Int, maximum(total_thickness) / dz))

    belt = zeros(Int, nx, ny, nz)

    for ix in 1:nx, iy in 1:ny
        zpos = 1

        for it in 1:nt
            dep = numeric_dep[:, ix, iy, it]
            thickness = sum(dep)

            thickness <= 0.0 && continue

            env = argmax(dep)
            nvox = max(1, round(Int, thickness / dz))
            ztop = min(nz, zpos + nvox - 1)

            belt[ix, iy, zpos:ztop] .= env

            zpos = ztop + 1
            zpos > nz && break
        end
    end

    return belt
end
# ~/~ end
# ~/~ begin <<docs/src/models/envcap.md#env-to-factory-prior>>[init]
function env_to_factory_prior_block(
        env_belt::Array{Int,3},
        mapping::Matrix{Float64})::Array{Float64,4}

    n_envs, n_factories = size(mapping)
    nx, ny, nz = size(env_belt)

    normalised = similar(mapping)

    for e in 1:n_envs
        row = max.(mapping[e, :], 0.0)
        s = sum(row)
        normalised[e, :] .= s > 0.0 ? row ./ s : fill(1.0 / n_factories, n_factories)
    end

    prior = zeros(Float64, n_factories, nx, ny, nz)

    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        e = env_belt[ix, iy, iz]

        if e == 0 || e > n_envs
            prior[:, ix, iy, iz] .= 1.0 / n_factories
        else
            prior[:, ix, iy, iz] .= normalised[e, :]
        end
    end

    return prior
end
# ~/~ end
# ~/~ begin <<docs/src/models/envcap.md#apply-prior-bias>>[init]
function apply_prior_bias!(
        ca::Matrix{Int},
        factory_prior::Array{Float64,4},
        sediment_height,
        depositional_resolution,
        ca_refinement::Float64,
        rng)

    ca_refinement == 0.0 && return

    _, nx, ny, nz = size(factory_prior)

    for ix in axes(ca, 1), iy in axes(ca, 2)
        f = ca[ix, iy]
        f == 0 && continue

        z = floor(Int, sediment_height[ix, iy] / depositional_resolution |> NoUnits) + 1
        z = clamp(z, 1, nz)

        p_prior = factory_prior[f, ix, iy, z]
        p_kill  = ca_refinement * (1.0 - p_prior)

        if rand(rng) < p_kill
            ca[ix, iy] = 0
        end
    end
end
# ~/~ end

end
# ~/~ end
