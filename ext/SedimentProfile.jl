module SedimentProfile

import CarboKitten.Visualization: sediment_profile, sediment_profile!, profile_plot!, coeval_lines!, profile_plot_coarse!, sediment_profile_generic!, facies_colormap
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, DataSlice
using CarboKitten.Skeleton: skeleton
using Makie
using GeometryBasics
using Unitful
using Statistics: mean
using Base: dropdims
export sediment_profile_generic!
const Time = typeof(1.0u"Myr")
const na = [CartesianIndex()]
const FACIES_COLORS = [:white, :dodgerblue, :hotpink, :lightyellow, :gold, :seagreen]



thickness_field(d::DataSlice) =
    hasproperty(d, :compacted_thickness) ? d.compacted_thickness : d.sediment_thickness

# -----------------------------------------------------------
# Safe aligned time
# -----------------------------------------------------------

aligned_time(header::Header, data::DataSlice) = begin
    thick = thickness_field(data)
    Nt1 = size(thick, 2)
    @assert Nt1 + 1 == length(header.axes.t)
    header.axes.t[1:Nt1]
end



# -----------------------------------------------------------
# Safe elevation
# -----------------------------------------------------------

function elevation(h::Header, d::DataSlice)

    thick = thickness_field(d)
    Nt1 = size(thick, 2)

    base = h.initial_topography[d.slice..., na]

    # interpretation-space view: final subsidence reference frame
    subs = d.cumulative_subsidence[:, Nt1]

    el = base .+ thick .- subs

    @assert size(el) == (length(h.axes.x), Nt1)

    return el
end



# -----------------------------------------------------------
# Mesh utilities
# -----------------------------------------------------------
"""
explode_quad_vertices(v)

Takes a three dimensional array representing a grid of vertices. This function duplicates these
vertices in the vertical direction, so that an amount of sediment can be given a single color.

Returns a tuple of vertices and faces (triangles), suitable for plotting with Makie's `mesh`
function.
"""

function explode_quad_vertices(v::Array{Float64,3})
    w, h, d = size(v)
    points = zeros(Float64, w, h - 1, 2, d)
    @views points[:,:,1,:] = v[:,1:end-1,:]
    @views points[:,:,2,:] = v[:,2:end,:]

    n_vertices = 2*w*(h-1)
    idx = reshape(1:n_vertices, w, h-1, 2)

    v1 = reshape(idx[1:end-1, :, 1], :)
    v2 = reshape(idx[2:end, :, 1], :)
    v3 = reshape(idx[2:end, :, 2], :)
    v4 = reshape(idx[1:end-1, :, 2], :)

    faces = vcat(hcat(v1,v2,v3), hcat(v1,v3,v4))

    return reshape(points, n_vertices, d), faces
end

# -----------------------------------------------------------
# Safe unconformities
# -----------------------------------------------------------

function plot_unconformities(ax::Axis, header::Header, data::DataSlice, ::Nothing; kwargs...) end

plot_unconformities(ax::Axis, header::Header, data::DataSlice, b::Bool; kwargs...) =
    b ? plot_unconformities(ax, header, data, 10; kwargs...) : nothing

function plot_unconformities(ax::Axis, header::Header, data::DataSlice, mw::Int; kwargs...)
    x = header.axes.x |> in_units_of(u"km")
    Nt1 = size(thickness_field(data), 2)


    ξ   = elevation(header, data)

sea  = header.sea_level[1:Nt1]

wd = sea' .- ξ

    @assert size(wd) == size(ξ)

    hi = skeleton(wd .> 0u"m", minwidth=mw)
    isempty(hi[1]) && return

    verts = [(x[i], ξ[i,j] |> in_units_of(u"m")) for (i,j) in hi[1]]
    segs = vec(permutedims(verts[hi[2]]))
    linesegments!(ax, segs;
    label = "Unconformities",
    kwargs...
)

end

# -----------------------------------------------------------
# Safe coeval lines
# -----------------------------------------------------------


function coeval_lines!(ax::Axis, header::Header, data::DataSlice, flag::Bool)
    flag && coeval_lines!(ax, header, data, (4,8))
end

function coeval_lines!(ax::Axis, header::Header, data::DataSlice, ts::Vector{Time}; kwargs...)
    t = aligned_time(header, data)
    idx = clamp.(searchsortedfirst.(Ref(t), ts), 1, length(t))
    coeval_lines!(ax, header, data, idx; kwargs...)
end

function coeval_lines!(ax::Axis, header::Header, data::DataSlice, ticset::Tuple{Int,Int})
    Nt1 = size(thickness_field(data), 2)
    major = round.(Int, range(1, Nt1, length=ticset[1]+1))
    minor = setdiff(round.(Int, range(1, Nt1, length=ticset[2]+1)), major)

    coeval_lines!(ax, header, data, minor;  color=:black, linewidth=1, linestyle=:dot)
    coeval_lines!(ax, header, data, major;  color=:black, linewidth=1, linestyle=:solid)
end

function coeval_lines!(ax::Axis, header::Header, data::DataSlice, idx::Vector{Int}; kwargs...)
    x = header.axes.x |> in_units_of(u"km")
    h = elevation(header, data) |> in_units_of(u"m")

    idx = clamp.(idx, 1, size(h,2))
    for t in idx
        lines!(ax, x, h[:,t];
    label = "Coeval lines",
    kwargs...
)

    end
end

# -----------------------------------------------------------
# Safe profile plot — no broadcast errors ever
# -----------------------------------------------------------


function profile_plot!(ax::Axis, header::Header, data::DataSlice; color, label = nothing, mesh_args...)

    x = header.axes.x |> in_units_of(u"km")
    n_f, n_x, n_t = size(data.production)
    ξ = elevation(header, data) |> in_units_of(u"m")

    verts = zeros(Float64, n_x, n_t+1, 2)
    verts[:,:,1] .= reshape(x, n_x, 1)
    verts[:,1:n_t,2] .= ξ
    verts[:,n_t+1,2] .= verts[:,n_t,2]

    v, f = explode_quad_vertices(verts)

    @assert size(color) == (n_x, n_t)
    c = repeat(vec(color), inner=2)

    return mesh!(ax, v, f; color = c, label = label, mesh_args...)


end

profile_plot!(f::F, ax::Axis, header::Header, data::DataSlice; mesh_args...) where F =
    profile_plot!(ax, header, data; color=f.(eachslice(data.deposition, dims=(2,3))), mesh_args...)

# -----------------------------------------------------------
# MAIN PROFILE
# -----------------------------------------------------------
"""
sediment_profile!(ax, header, data; show_unconformities)

Plot the sediment profile, choosing colour by dominant facies type (argmax). Unconformaties
are shown when the sediment is subaerially exposed (even if sediment is still deposited
due to a set intertidal zone).
"""
function sediment_profile!(ax::Axis, header::Header, data::DataSlice;
    show_unconformities=true,
    show_coeval_lines=true,
    show_sealevel=true,
    facies_colors=nothing
)
    n_f = size(data.production,1)
    dep = data.deposition

    dom = argmax.(eachslice(dep, dims=(2,3)))
    tot = dropdims(sum(dep; dims=1), dims=1)

    domf = Float32.(dom)
    domf[tot .<= 0u"m"] .= NaN32

    cmap = facies_colormap(n_f; facies_colors=facies_colors, include_nodeposit=false)

    hm = profile_plot!(ax, header, data;
    color = domf,
    alpha = 1.0,
    colormap = cmap,
    colorrange=(0.5f0, n_f+0.5f0),
    label = "Dominant facies"
)


    # -------------------------------
    # Fixed coeval-lines dispatch
    # -------------------------------
    if show_coeval_lines isa Bool
        show_coeval_lines && coeval_lines!(ax, header, data, (4,8))

    elseif show_coeval_lines isa Tuple{Int,Int}
        coeval_lines!(ax, header, data, show_coeval_lines)

    elseif show_coeval_lines isa Vector{Time}
        coeval_lines!(ax, header, data, show_coeval_lines)

    elseif show_coeval_lines isa Vector{Int}
        coeval_lines!(ax, header, data, show_coeval_lines)

    else
        @warn "Unknown coeval_lines type $(typeof(show_coeval_lines)) — skipping."
    end

    # -------------------------------
    # Unconformities
    # -------------------------------
    show_unconformities !== false &&
        plot_unconformities(ax, header, data, show_unconformities;
            color=:white, linestyle=:dash)

    # -------------------------------
    # Sea level
    # -------------------------------
    if show_sealevel
        sl = header.sea_level[end] |> in_units_of(u"m")
        hlines!(ax, sl;
    color = :lightblue,
    linewidth = 5,
    label = "Sea level"
)

    end

    ax.title = "sediment profile"
    return ax, hm
end

function sediment_profile(header::Header, data::DataSlice; kwargs...)
    fig = Figure(size=(1000,600))
    ax  = Axis(fig[1,1])
    ax, hm = sediment_profile!(ax, header, data; kwargs...)
    return fig, ax, hm
end



function profile_plot_coarse!(
    ax::Axis,
    header::Header,
    data::DataSlice;
    color,
    Δi::Int,
    label = nothing,
    coarse_colors = nothing,
    mesh_args...
)

    # coarse x axis
    x_fine = header.axes.x |> in_units_of(u"km")
    Nx_fine = length(x_fine)

    Nx_coarse = size(color, 1)
    n_t       = size(color, 2)

    @assert Nx_coarse == ceil(Int, Nx_fine / Δi)

    # block centers
    x_coarse = [
        mean(x_fine[(i-1)*Δi+1:min(i*Δi, Nx_fine)])
        for i in 1:Nx_coarse
    ]

    ξ = elevation(header, data) |> in_units_of(u"m")
    ξ_coarse = [
        mean(ξ[(i-1)*Δi+1:min(i*Δi, Nx_fine), :]; dims=1)
        for i in 1:Nx_coarse
    ] |> x -> reduce(vcat, x)

    verts = zeros(Float64, Nx_coarse, n_t+1, 2)
    verts[:,:,1] .= reshape(x_coarse, Nx_coarse, 1)
    verts[:,1:n_t,2] .= ξ_coarse
    verts[:,n_t+1,2] .= verts[:,n_t,2]

    v, f = explode_quad_vertices(verts)

    c = repeat(vec(color), inner=2)

    n_f = size(data.production, 1)

cols = coarse_colors === nothing ?
    Makie.wong_colors() :
    Makie.to_color.(coarse_colors)

cmap = cgrad(cols, n_f; categorical=true)

mesh!(ax, v, f;
    color = c,
    colormap = cmap,
    colorrange = (0.5f0, n_f + 0.5f0),
    label = label,
    mesh_args...
)
end

function elevation_physical(h::Header, d::DataSlice)

    thick = thickness_field(d)
    Nt1 = size(thick, 2)

    base = h.initial_topography[d.slice..., na]

    # process-space view: time-dependent subsidence
    subs = d.cumulative_subsidence[:, 1:Nt1]

    el = base .+ thick .- subs

    return el
end


function sediment_profile_generic!(
    ax::Axis,
    header::Header,
    data::DataSlice;
    selector::Function,
    category_names::Vector{String},
    category_colors,
    facies_colors = nothing,
    nondep_eps = 1e-12,
    show_unconformities = true,
    show_coeval_lines = true,
    show_sealevel = true,
    colorbar_label = ""
)

    dep = data.deposition
    n_x = size(dep, 2)
    n_t = size(dep, 3)

    # Precompute elevation once (units preserved)
    ξ_phys = elevation_physical(header, data)

    # Use matching sea-level segment
    sea = header.sea_level[1:n_t]

    dom = zeros(Float32, n_x, n_t)

    for ix in 1:n_x, it in 1:n_t
        strength = dep[:, ix, it]
        total = sum(ustrip.(strength))

        if total ≤ nondep_eps
            dom[ix, it] = 0.0f0
        else
            # water depth (positive when submerged)
            wd = sea[it] - ξ_phys[ix, it]
			wd_m = ustrip(wd |> in_units_of(u"m"))

            e_kwpm = begin
    if hasproperty(data, :wave_energy)
        Float64(data.wave_energy[ix, it])
    else
        NaN
    end
end

		dom[ix, it] = try
			Float32(selector(strength, wd_m, e_kwpm))
		catch err
			err isa MethodError || rethrow()
			Float32(selector(strength, wd_m))  # backward compatible
		end
        end
    end

    category_colors = facies_colors === nothing ? category_colors : facies_colors
    cmap = facies_colormap(length(category_names); facies_colors=category_colors, include_nodeposit=false)

    hm = profile_plot!(ax, header, data;
        color = dom,
        colormap = cmap,
        colorrange = (-0.5f0, length(category_names) - 0.5f0),
        label = colorbar_label
    )

    if show_coeval_lines isa Bool
        show_coeval_lines && coeval_lines!(ax, header, data, (4,8))
    elseif show_coeval_lines isa Tuple{Int,Int}
        coeval_lines!(ax, header, data, show_coeval_lines)
    elseif show_coeval_lines isa Vector
        coeval_lines!(ax, header, data, show_coeval_lines)
    end

    show_unconformities !== false &&
        plot_unconformities(ax, header, data, show_unconformities;
            color = :white, linestyle = :dash)

    if show_sealevel
        sl = header.sea_level[end] |> in_units_of(u"m")
        hlines!(ax, sl; color = :lightblue, linewidth = 5)
    end

    return ax, hm
end

end # module
