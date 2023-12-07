module climateconfig
export climate,model1

struct climate
    P::Float64
    T::Float64
    pco2::Float64
end

model1 = [climate(1000.0,288.0,10^(-1.5))]
end