# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/fig8-sealevel.jl>>[init]
#| creates: data/bs92-sealevel-curve.csv
#| requires: data/bs92-sealevel-input.png

module Script
    using Images
    using DataFrames
    using CSV

    function main()
        img = load("data/bs92-sealevel-input.png")
        img_gray = Gray.(img)
        signal = 1.0 .- channelview(img_gray)
        signal ./= sum(signal; dims=[1])
        (n_y, n_x) = size(signal)
        y = sum(signal .* (1:n_y); dims=[1]) / n_y * 200.0
        df = DataFrame(
            time = LinRange(0.0, 80_000.0, n_x),
            depth = y[1, :])
        CSV.write("data/bs92-sealevel-curve.csv", df)
    end
end

Script.main()
# ~/~ end