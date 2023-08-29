---
title: CarboKitten composite model
subtitle: CA + production + transport
---

``` {.julia #carbocat-composite}
function run()
    # Run the CA for 10 generations as a warm-up
    species = Iterators.drop(CA.run(Reflected{2}, rand(0:3, 50, 50), 3), 10)
    # Initial depth runs from 0 to 150
    height = repeat(collect(0:49) .* 3.0, 1, 50)
    sealevel(t) = 0.0
    Δt = 1000.0

    for (time_index, gen) in enumerate(species)
        t = time_index * Δt
        production = Δt .* Production.production_rate(2000.0, Config.model1, height .- sealevel(t))

    end
end
```
