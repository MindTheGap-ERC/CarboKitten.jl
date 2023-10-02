# CarboKitten
**Modeling Carbonate Platforms in Julia**

[![Entangled badge](https://img.shields.io/badge/entangled-Use%20the%20source!-%2300aeff)](https://entangled.github.io/)

## Project layout


## Running

Start the Julia REPL, and get into Pkg mode by pressing `]`. You may activate the package environment using `activate .` and then install the dependencies using `instantiate`. These steps only need to be run once.

```julia
pkg> activate .
pkg> instantiate
```

If you want to start a REPL with the correct environment already activated, use the `--project=.` flag. Use the `-t` flag to enable processing in multiple threads.

```shell
julia --project=. -t 4
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
include("examples/bosscher-schlager-1992.jl")
```

After that, you may edit an example and rerun.

## Development
While developing, you'll need to run the Entangled watch daemon to keep Markdown and Julia code synchronized.

```shell
entangled watch
```

The documentation is generated using `Documenter.jl`. The most efficient way to serve this documentation and have it update upon changes, is to run `LiveServer` from the Julia REPL or

```shell
make serve-docs
```

To generate the more expensive figures (actually resulting from simulation etc.), you need to run a `DaemonMode` process in the background.

```shell
make run-daemon
```

After that, you can run

```shell
make figures
```

To run simulations and plot figures depending on those.

## License

Copyright University of Utrecht and Netherlands eScience Center 2023. This repository is licensed under the Apache v2 license, see LICENSE.

## References

Code in this repository is based on

- Burgess 2013
- Bosscher and Schlager 1992
