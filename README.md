# CarboKitten

**Modeling Carbonate Platforms in Julia**

## Running

Start the Julia REPL, and get into Pkg mode by pressing `]`. You may activate the package environment using `activate .` and then install the dependencies using `instantiate`. These steps only need to be run once.

```julia
pkg> activate .
pkg> instantiate
```

If you want to start a REPL with the correct environment already activated, use the `--project=.` flag.

```shell
julia --project=.
```

### Examples

You'll get the best experience by running examples from the Julia REPL. There are however also some example scripts, that should work stand-alone.

```shell
julia --project=. examples/bosscher-schlager-1992.jl
```

However, it is more efficient to run them from the REPL. Either run,

```shell
julia --project=.
```

or start the REPL from VS Code. In the REPL you can run

```julia
include("examples/bosscher-sclager-1992.jl")
```

After that, you may edit an example and rerun.

## License

Copyright University of Utrecht and Netherlands eScience Center 2023. This repository is licensed under the Apache v2 license, see LICENSE.

## References

Code in this repository is based on

- Burgess 2013
- Bosscher and Schlager 1992
