# CarboKitten
**Modeling Carbonate Platforms in Julia**

[![Entangled badge](https://img.shields.io/badge/entangled-Use%20the%20source!-%2300aeff)](https://entangled.github.io/)

## Project layout

```
.
├── data                # data files
├── docs                # documentation
│   ├── make.jl         # docs build script
│   ├── Manifest.toml   # 
│   ├── Project.toml    # dependencies for building docs
│   └── src             # markdown source for docs
├── entangled.toml      # entangled config
├── examples            # example scripts
├── Makefile            # command-line short hands
├── Manifest.toml       #
├── Project.toml        # project dependencies
├── README.md           # 
└── src                 # tangled library source
```

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

You'll get the best experience by running examples from the Julia REPL. There are however also some example scripts that should work stand-alone.

```shell
julia --project=. examples/ca-with-prod.jl
```
This command will write the output in the HDF5 format into the `data` folder. You can check that output is written there after executing this command.

However, it is more efficient to run them from the REPL. Either run,

```shell
julia --project=.
```

or start the REPL from VS Code. In the REPL you can run

```julia
include("examples/ca-with-prod.jl")
```

After that, you may edit an example and rerun.

## Development

### Global dependencies

CarboKitten has some dependencies that are only needed for developing, but not for using the library. You may want to install those in the global environment:

```shell
julia -e 'using Pkg; Pkg.add.(["DaemonMode", "LiveServer", "Revise"])'
```

### Entangled

While developing, you'll need to run the [Entangled](https://entangled.github.io/) watch daemon to keep Markdown and Julia code synchronized.

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

## Documentation

To generate the documentation, run `julia`.

```
pkg> activate docs
pkg> instantiate
julia> using LiveServer
julia> servedocs()
```

The figures are generated seperately.

## References

Code in this repository is based on

- Burgess, P. M. (2013). [CarboCAT: A cellular automata model of heterogeneous carbonate strata](https://www.sciencedirect.com/science/article/pii/S0098300411002949). Computers & geosciences, 53, 129-140.
- Bosscher, H., & Schlager, W. (1992). [Computer simulation of reef growth](https://doi.org/10.1111/j.1365-3091.1992.tb02130.x). Sedimentology, 39(3), 503-512.

## Authors

Lead engineer: __Johan Hidding__  
The Netherlands eScience Center  
email: j.hidding [at] esciencecenter.nl   
Web page: [www.esciencecenter.nl/team/johan-hidding-msc/](https://www.esciencecenter.nl/team/johan-hidding-msc/)  
ORCID: [0000-0002-7550-1796](https://orcid.org/0000-0002-7550-1796)

Original author: __Peter Burgess__  
University of Liverpool  
Web page: [www.liverpool.ac.uk/environmental-sciences/staff/peter-burgess](https://www.liverpool.ac.uk/environmental-sciences/staff/peter-burgess/)

Project lead: __Emilia Jarochowska__  
Utrecht University  
email: e.b.jarochowska [at] uu.nl  
Web page: [www.uu.nl/staff/EBJarochowska](https://www.uu.nl/staff/EBJarochowska)  
ORCID: [0000-0001-8937-9405](https://orcid.org/0000-0001-8937-9405)

**Other team members:**

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Xianyi Liu__  
Utrecht University  
email: x.liu6 [at] uu.nl  
Web page: [www.uu.nl/staff/XLiu6](https://www.uu.nl/staff/XLiu6)  
ORCID: 

__Hanno Spreeuw__  
The Netherlands eScience Center  
email: h.spreeuw [at] esciencecenter.nl  
Web page: [www.esciencecenter.nl/team/dr-hanno-spreeuw/](https://www.esciencecenter.nl/team/dr-hanno-spreeuw)  
ORCID: [0000-0002-5057-0322](https://orcid.org/0000-0002-5057-0322)

__David De Vleeschouwer__  
Westfälische Wilhelms-Universität Münster  
Web page: [www.uni-muenster.de/GeoPalaeontologie/erdsystemforschung/staff/DeVleeschouwer](https://www.uni-muenster.de/GeoPalaeontologie/erdsystemforschung/staff/DeVleeschouwer.html)  
ORCID: [0000-0002-3323-807X](https://orcid.org/0000-0002-3323-807X)

## Copyright

Copyright 2023-2024 The Netherlands eScience Center and Utrecht University

## License

Apache 2.0 License, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.

