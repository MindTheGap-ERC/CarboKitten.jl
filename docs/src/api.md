# API Documentation

```@docs
run_model
box_axes
CarboKitten.Box
TimeProperties
time_axis
n_steps
```

## Production

```@autodocs
Modules = [CarboKitten.Production]
```

## Transport

```@autodocs
Modules = [CarboKitten.Transport.Advection]
```

## Export

```@autodocs
Modules = [CarboKitten.Export]
```

## Components

```@autodocs
Modules = [CarboKitten.Components.CellularAutomaton, CarboKitten.Components.ActiveLayer]
```

## Output

```@autodocs
Modules = [CarboKitten.Output.Abstract, CarboKitten.Output.H5Writer, CarboKitten.Output.MemoryWriter]
```

## Utility

```@autodocs
Modules = [CarboKitten.Utility, CarboKitten.Stencil]
```

## Submodules

```@autodocs
Modules = [
  CarboKitten.Denudation, CarboKitten.Denudation.Abstract,
  CarboKitten.Denudation.NoDenudationMod
  ]
```

## Algorithms

```@autodocs
Modules = [
  CarboKitten.Algorithms.EnumerateSeq,
  CarboKitten.Algorithms.Skeleton,
  CarboKitten.Algorithms.RangeFinder,
  CarboKitten.Algorithms.StratigraphicColumn,
  CarboKitten.Algorithms.Romberg
]
```
