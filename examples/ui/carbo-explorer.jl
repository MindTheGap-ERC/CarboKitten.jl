module Explorer

using GLMakie
using Unitful

using CarboKitten.Export: Header, DataSlice, read_slice, age_depth_model
using CarboKitten.Visualization: sediment_profile!, production_curve!, wheeler_diagram!
using CarboKitten.Model.ALCAPS: Input

const INPUT = Input()
const Amount = typeof(1.0u"m")

struct DataColumn
    slice::NTuple{2,Int} 
    disintegration::Array{Amount,2}
    production::Array{Amount,2}
    deposition::Array{Amount,2}
    sediment_elevation::Array{Amount,1}
end

function stratigraphic_column(header::Header, data::DataColumn, facies::Int)
    n_times = length(header.axes.t) - 1
    sc = zeros(typeof(1.0u"m"), n_times)

    for ts = 1:n_times
        acc = data.deposition[facies, ts] - data.disintegration[facies, ts]
        if acc > 0.0u"m"
            sc[ts]= acc
            continue
        end
        ts_down = ts - 1
        while acc < 0.0u"m"
            ts_down < 1 && break
            if -acc < sc[ts_down]
                sc[ts_down] -= acc
                break
            else
                acc += sc[ts_down]
                sc[ts_down] = 0.0u"m"
            end
            ts_down -= 1
        end
    end

    sc
end

function scdata(header::Header, data::DataColumn)
    n_facies = size(data.production)[1]
    n_times = length(header.axes.t) - 1
    sc = zeros(Float64, n_facies, n_times)
    for f = 1:n_facies
        sc[f,:] = stratigraphic_column(header, data, f) / u"m"
    end
    
    colormax(d) = getindex.(argmax(d; dims=1)[1, :], 1)
    adm = age_depth_model(data.sediment_elevation)

    return (ys_low=adm[1:end-1]/u"m", ys_high=adm[2:end]/u"m", color=Makie.wong_colors()[colormax(sc)[1:end-1]])
end

function main(filename)
    y = 25
    header, data = read_slice(filename, :, y)

    fig = Figure(size=(1000, 1080))
    ax_sp = Axis(fig[1,1:2])
    sl_x = Slider(fig[2,1:2], range=1:length(header.axes.x), startvalue=50)
    ax_wd1 = Axis(fig[4,1])
    ax_wd2 = Axis(fig[4,2])
    ax_pc = Axis(fig[1,3])

    column = lift(sl_x.value) do x
        DataColumn((x, y), data.disintegration[:,x,:],
                   data.production[:,x,:],
                   data.deposition[:,x,:],
                   data.sediment_elevation[x,:])
    end

    _scdata = lift(c->scdata(header, c), column)
    _sc_low = lift(c->c.ys_low, _scdata)
    _sc_hig = lift(c->c.ys_high, _scdata)
    _sc_col = lift(c->c.color, _scdata)

    ax_sc = Axis(fig[4,4], width=50)
    hspan!(ax_sc, _sc_low, _sc_hig; color=_sc_col)
    _adm = lift(column) do c
        age_depth_model(c.sediment_elevation) / u"m"
    end

    ax_adm = Axis(fig[4,3])
    lines!(ax_adm, header.axes.t[1:end-1] / u"Myr", _adm)

    sediment_profile!(ax_sp, header, data)
    production_curve!(ax_pc, INPUT)
    
    (wd1, wd2) = wheeler_diagram!(ax_wd1, ax_wd2, header, data)
    for ax in [ax_wd1, ax_wd2, ax_sp]
        vlines!(ax, lift(x -> header.axes.x[x]/u"km", sl_x.value); linewidth=3, color=(:white, 0.5))
        vlines!(ax, lift(x -> header.axes.x[x]/u"km", sl_x.value); linewidth=1, color=:black)
    end
    
	Colorbar(fig[3,1], wd1; vertical=false, label="sediment accumulation [m/Myr]")
	Colorbar(fig[3,2], wd2; vertical=false, ticks=1:3, label="dominant facies")

    fig
end

end

Explorer.main("data/alcaps2.h5")

