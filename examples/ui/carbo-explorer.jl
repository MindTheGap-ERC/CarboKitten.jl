module Explorer

using GLMakie
using Unitful
using HDF5

using CarboKitten.Export: Header, DataSlice, DataColumn, read_slice, read_header, age_depth_model
using CarboKitten.Visualization: sediment_profile!, production_curve!, wheeler_diagram!, stratigraphic_column!

const Amount = typeof(1.0u"m")

function main(filename)
    fig = Figure(size=(1000, 1080))
    h5open(filename, "r") do fid
        header = read_header(fid)
        y = length(header.axes.y) ÷ 2 + 1
        data = read_slice(fid, :, y)

        ax_sp = Axis(fig[1, 1:2])
        sl_x = Slider(fig[2, 1:2], range=1:length(header.axes.x), startvalue=50)
        ax_wd1 = Axis(fig[4, 1])
        ax_wd2 = Axis(fig[4, 2])
        ax_pc = Axis(fig[1, 3])

        column = lift(sl_x.value) do x
            DataColumn((x, y), data.disintegration[:, x, :],
                data.production[:, x, :],
                data.deposition[:, x, :],
                data.sediment_elevation[x, :])
        end

        ax_sc = Axis(fig[4, 4]; width=100)
        stratigraphic_column!(ax_sc, header, column)

        _adm = lift(column) do c
            age_depth_model(c.sediment_elevation) / u"m"
        end

        ax_adm = Axis(fig[4, 3])
        lines!(ax_adm, header.axes.t[1:end] / u"Myr", _adm)

        sediment_profile!(ax_sp, header, data)
        production_curve!(ax_pc, fid["input"])

        (wd1, wd2) = wheeler_diagram!(ax_wd1, ax_wd2, header, data; smooth_size=(1,1))
        for ax in [ax_wd1, ax_wd2, ax_sp]
            vlines!(ax, lift(x -> header.axes.x[x] / u"km", sl_x.value); linewidth=3, color=(:white, 0.5))
            vlines!(ax, lift(x -> header.axes.x[x] / u"km", sl_x.value); linewidth=1, color=:black)
        end

        Colorbar(fig[3, 1], wd1; vertical=false, label="sediment accumulation [m/Myr]")
        Colorbar(fig[3, 2], wd2; vertical=false, ticks=1:3, label="dominant facies")
    end
    fig
end

end

Explorer.main("data/output/alcap2.h5")

