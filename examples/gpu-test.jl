module GPUTest

using oneAPI
using KernelAbstractions

using Unitful
using CarboKitten
using CarboKitten: AbstractInput
using CarboKitten.BoundaryTrait: Boundary, Periodic, Reflected, Coast
using CarboKitten.Models: ALCAP

using StaticArrays

@inline modflip(a, l) =
    let b = mod1(a, 2l)
        b > l ? 2l - b + 1 : b
    end

@inline get_bounded(::Type{Periodic{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[mod1.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Reflected{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[modflip.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Coast}, a, i) =
    checkbounds(Bool, a, i) ? a[i] : a[modflip(i[1], size(a)[1]), mod1(i[2], size(a)[2])]

@kernel function life(::Type{BT}, a, b) where {BT}
    i = @index(Global, Cartesian)

    n = 0
    for j in CartesianIndices((3, 3))
        n += get_bounded(BT, a, i + j - CartesianIndex(1, 1))
    end
    n -= a[i]

    b[i] = n == 3 || n == 2 && a[i] == 1
end

@kernel function ca2d5sq(::Type{BT}, in, out, precedence) where {BT}
    i = @index(Global, Cartesian)

    for f in precedence
        n = 0
        for j in CartesianIndices((5, 5))
            n += (get_bounded(BT, in, i + j - CartesianIndex(3, 3)) == f)
        end

        if in[i] == f
            n -= 1
            out[i] = 4 <= n && n <= 10 ? f : 0
            break
        else
            if 6 <= n && n <= 10
                out[i] = f
                break
            end
        end
    end
end

@inline function production_rate(insolation, maximum_growth_rate, saturation_intensity, extinction_coefficient, water_depth)
    gₘ = maximum_growth_rate
    I = insolation / saturation_intensity
    x = water_depth * extinction_coefficient
    return x > 0.0f0 ? gₘ * tanh(I * exp(-x)) : 0.0f0
end

function production(input::AbstractInput)
    :(@kernel function (water_depth, facies_map, production_out)
        I_0 = $(input.insolation |> in_units_of(u"W/m^2") |> Float32)
        dt = $(input.time.Δt |> in_units_of(u"Myr") |> Float32)
        g_m = $([f.maximum_growth_rate |> in_units_of(u"m/Myr") |> Float32 for f in input.facies])
        I_k = $([f.saturation_intensity |> in_units_of(u"W/m^2") |> Float32 for f in input.facies])
        k = $([f.extinction_coefficient |> in_units_of(u"1/m") |> Float32 for f in input.facies])
        i = @index(Global, Cartesian)

        f = facies_map[i]
        if f != 0
            wd = water_depth[i]
            pr = production_rate(I_0, g_m[f], I_k[f], k[f], wd) * dt
            production_out[f, i[1], i[2]] = min(pr, max(0.0f0, wd))
        end
    end)
end

function main()
    dev = CPU()
    facies = [
        ALCAP.Facies(
            maximum_growth_rate=500u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=50.0u"m/yr"),
        ALCAP.Facies(
            maximum_growth_rate=400u"m/Myr",
            extinction_coefficient=0.1u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=25.0u"m/yr"),
        ALCAP.Facies(
            maximum_growth_rate=100u"m/Myr",
            extinction_coefficient=0.005u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=12.5u"m/yr")
    ]

    input = ALCAP.Input(
        box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
        time=TimeProperties(
            Δt=0.0002u"Myr",
            steps=5000),
        output=Dict(
            :topography => OutputSpec(slice=(:,:), write_interval=10),
            :profile => OutputSpec(slice=(:, 25), write_interval=1)),
        ca_interval=1,
        initial_topography=(x, y) -> -x / 300.0,
        sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
        subsidence_rate=50.0u"m/Myr",
        disintegration_rate=50.0u"m/Myr",
        insolation=400.0u"W/m^2",
        sediment_buffer_size=50,
        depositional_resolution=0.5u"m",
        facies = facies
    )

    # a = oneArray(rand(Int32(0):Int32(3), 256, 256))
    # b = oneArray(zeros(Int32, 256, 256))
    # total_sediment = oneArray(zeros(Float32, 3, 256, 256))
    # wd = oneArray(fill(10.0f0, 256, 256))
    a = rand(Int32(0):Int32(3), 256, 256)
    b = zeros(Int32, 256, 256)
    total_sediment = zeros(Float32, 3, 256, 256)
    wd = fill(10.0f0, 256, 256)

    precedence = collect(1:3)
    p_expr = production(input)
    # print(p_expr)
    production_k = eval(p_expr)

    backend = get_backend(a)
    for _ in 1:1000
        ca2d5sq(backend, 256)(Periodic{2}, a, b, precedence, ndrange=size(a))
        KernelAbstractions.synchronize(backend)
        circshift!(precedence, 1)
        a, b = b, a

        production_k(backend, 256)(wd, a, total_sediment, ndrange=size(a))

        KernelAbstractions.synchronize(backend)
    end
    Array(total_sediment)
end

end
