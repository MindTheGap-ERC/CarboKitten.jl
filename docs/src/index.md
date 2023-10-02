# CarboKitten
**Modeling Carbonate Platforms in Julia**

[![Entangled badge](https://img.shields.io/badge/entangled-Use%20the%20source!-%2300aeff)](https://entangled.github.io/)

![](fig/b13-crosssection.png)


```@contents
```

## Julia Quickstarter

This code is written in [Julia](https://julia-lang.org/). You may want to check out the following references:

- [Julia Documentation](https://docs.julialang.org/en/v1/)
- [Tutorial on Julia for Science and Engineering](https://www.matecdev.com/posts/julia-tutorial-science-engineering.html)

There are several ways to work with Julia that may be a bit different from what you're used to, if that is Matlab, Python or R.

### REPL
The most basic way to work in Julia, is to start the REPL (Read Eval Print Loop).

```shell
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```

From here you may use CarboKitten `using CarboKitten` and run any of the code inside. To work with CarboKitten efficiently, you may want to load `Revise`. Revise auto-detects changes to loaded code and makes it easy to rerun.

Additionally you should learn how to work with Julia packages. If you want to experiment with things, try to create a new environment in an empty folder and add CarboKitten as a `dev` dependency:

```
pkg> dev <path to CarboKitten>
```

### Jupyter
You can run Julia code from Jupyter if you install the Julia kernel. Press `]` in the REPL to get into Pkg-mode, the prompt will change

```shell
(CarboKitten) pkg>
```

You may install the `IJulia` kernel with `add IJulia`.

### Pluto
An alternative notebook interface is called `Pluto`.

- Pluto is **reactive**: changes to code cells automatically update downstream dependencies.
- Pluto notebooks are written to regular Julia files and can be run independent from Pluto.
- The user interface of Pluto is slightly less mature than Jupyter

In Pkg-mode say `add Pluto`.

```shell
julia> using Pluto

julia> Pluto.run()
[ Info: Loading...
┌ Info: 
└ Opening http://localhost:1234/?secret=xyzxyzzy in your default browser... ~ have fun!
┌ Info: 
│ Press Ctrl+C in this terminal to stop Pluto
└ 
```

