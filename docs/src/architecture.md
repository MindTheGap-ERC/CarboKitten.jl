# Architecture

CarboKitten is modular in design. This gives us the advantage of few repetitions in the code, but it also means a little study getting in to developing CarboKitten.

## Models

A model is a module that must have certain data structures and methods implemented.

- `struct Input <: AbstractInput`: contains all input parameters
- `struct State <: AbstractState`: contains all run-time state
- `function initial_state(input)`: creates an initial `State` from an `Input`
- `function step!(input)(state)`: curried function to advance state and return a `H5Writer.DataFrame`
- `function write_header(fid, input)`: writes meta-data to an open HDF5 file

For most models the `Input` and `State` structures are generated from a hierarchy of components, which is explained below.

### `step!`

The `step!` method deserves some extra attention: this is where the high-level logic of a model goes. Running a model is nothing but a repeated application of `step!(input)` on the state. Note that the `step!` function is [curried](https://en.wikipedia.org/wiki/Currying): the `input` and `state` parameters are given on separate occasions. This allows the model to prepare later execution on `state` in a more efficient manner. For example, we may allocate memory or prepare quantities derived from input variables.

## Components

A model in CarboKitten is composed of components. Each component extends the `Input` and/or `State` data structures with members and can also provide functions that work on those elements. Components can inherit from each other. For instance, the `WaterDepth` component inherits from `Boxes` component to have a notion of the grid size, and the `TimeIntegration` component to obtain a time coordinate from the state.

A component has the following syntax:

```julia
@compose module MyComponent
    @mixin DependencyA, DependencyB
    using ..Common

    @kwdef struct Input <: AbstractInput
        ...
    end

    ...
end
```

Here the `Common` module is a normal Julia module living in the same scope as the component. This `Common` module contains common definitions used throughout the CarboKitten code base.

## H5Writer

The `H5Writer` component is special as it serves to execute a model. For example, if we have a model `MyModel` with all the required definitions,

```julia
@component module MyModel
    @mixin Tag, H5Writer, ...
    ...

    function step!(input)
        ...
        return function (state)
            ...
            return H5Writer.DataFrame(...)
        end
    end

    ...
end
```

We can run that model as follows:

```julia
H5Writer.run(Model{MyModel}, MyModel.Input(), "output.h5")
```

This will call `MyModel.initial_state` first, then repeatedly `MyModel.step!`, which should return a `H5Writer.DataFrame` for writing to HDF5. Here `Model` is a [value type](https://docs.julialang.org/en/v1/manual/types/#%22Value-types%22) (defined in `CarboKitten.Components.Common`).
