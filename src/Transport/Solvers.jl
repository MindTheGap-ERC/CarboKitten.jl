# ~/~ begin <<docs/src/finite-difference-transport.md#src/Transport/Solvers.jl>>[init]
module Solvers

using Unitful

function runge_kutta_4(box)
    k1 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k2 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k3 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k4 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    function (df, y, t, dt)
        k1 .= df(y, t)
        k2 .= df(y .+ dt/2 .* k1, t + dt/2)
        k3 .= df(y .+ dt/2 .* k2, t + dt/2)
        k4 .= df(y .+ dt .* k3, t + dt)
        y .+= (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) .* (dt/6)
    end
end

function forward_euler(df, y, t, dt)
    y .+= dt .* df(y, t)
end

end
# ~/~ end