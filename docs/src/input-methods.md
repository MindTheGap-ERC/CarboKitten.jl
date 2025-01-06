# Input Methods

Every model in CarboKitten has its own specific needs in terms of input. However, there are many common idioms used here.

## Sea level functions
Principally, the sea level curve is probably the input parameter that is most often of interest.

```julia
const INPUT = Input(
    ...
    sea_level = # what do I enter here?
    ...)
```

There are several ways that Julia allows us to enter functions.

### Inline functions
If the function is very simple, you can enter it as an inline anonymous function (also known as a lambda). The following would generate a sinusoid with an amplitude of 10m and a period of 100,000 years.

```julia
const INPUT = Input(
    ...
    sea_level = t -> 10.0u"m" * sin(2pi / 100.0u"kyr")
    ...)
```

### Stochastic functions
Functions can capture pre-generated data, so you can generate a stochastic sea-level curve for use in CarboKitten.

### Using a prepared dataset
CarboKitten ships with the Cenozoic sea level dataset from Miller 2020.

### From external tables
Using other packages like `DelimitedFiles`, `XLSX.jl` and `Interpolations.jl` you can import and interpolate your own data for use with CarboKitten.
