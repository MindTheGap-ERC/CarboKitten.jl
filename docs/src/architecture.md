# Architecture

CarboKitten is modular in design. This gives us the advantage of few repetitions in the code, but it also means a little study getting in to developing CarboKitten.

## Models

A model is a module that must have certain data structures and methods implemented.

- `struct Input <: AbstractInput`: contains all input parameters
- `struct State <: AbstractState`: contains all run-time state
- `function initial_state(input)`: creates an initial `State` from an `Input`
- `function step!(input)(state)`: curried function to advance state and return a `Frame`
- `function write_header(fid, input)`: writes meta-data to an open HDF5 file

For most models the `Input` and `State` structures are generated from a hierarchy of components, which is explained below.

### `step!`

The `step!` method deserves some extra attention: this is where the high-level logic of a model goes. Running a model is nothing but a repeated application of `step!(input)` on the state. Note that the `step!` function is [curried](https://en.wikipedia.org/wiki/Currying): the `input` and `state` parameters are given on separate occasions. This allows the model to prepare later execution on `state` in a more efficient manner. For example, we may allocate memory or prepare quantities derived from input variables.

The `step!` function should modify the state and return a `Frame`, which is defined as follows:

```julia
@kwdef struct Frame
    disintegration::Union{Array{Sediment,3},Nothing} = nothing   # facies, x, y
    production::Union{Array{Sediment,3},Nothing} = nothing
    deposition::Union{Array{Sediment,3},Nothing} = nothing
end
```

Here `disintegration` is the amount of sediment that was disintegrated, `production` the amount of sediment produced in situ, and `deposition` the amount of sediment deposited from the active layer. The only reason why this `Frame` is passed on, is for later storage. So the only routine that needs to understand the `Frame` structure is the `write_frame` method belonging to your output writer (currently `H5Writer` or `MemoryOutput`).

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

    @kwdef struct State <: AbstractState
        ...
    end

    ...
end
```

Here the `Common` module is a normal Julia module living in the same scope as the component. This `Common` module contains common definitions used throughout the CarboKitten code base.

### Module Mixins

CarboKitten uses a non-idiomatic method of abstraction to structure the components. The `@compose` macro is rather special. It takes care of a form of inheritance that is normally absent in the Julia language. All the modules that are mentioned in a `@mixin` clause will be merged into a single module, where all the structs that were defined in those modules are combined. This allows us to layer our code, separating concerns.

The problem is that this abstraction is neither idiomatic for those used to Julia nor those used to traditional object oriented languages. These layers together build two major important structures: `Input` and `State`. This includes any nested data structures (mainly `Facies` for the input data).

| Component  | Input                | Facies            | State             |
| ---------- | -------------------- | ----------------- | ----------------- |
| Boxes      | `box`                |                   |                   |
| FaciesBase | `facies`             | `label = Nothing` |                   |
| Time       | `time`               |                   | `step`            |
| WaterDepth | `sea_level`          |                   | `sediment_height` |
|            | `initial_topography` |                   |                   |
|            | `subsidence_rate`    |                   |                   |
| etc.       |                      |                   |                   |

The Julia language allows us to work with abstract interfaces by using dispatch, but it has no built-in method to compose fully formed concrete structures like this. Note that the `WaterDepth` component defines `subsidence_rate`, but to compute water depth from that, we need to know the time. This is why the `WaterDepth` component depends/inherits/mixins the `Time` component.

Each component is automatically documented with a graph that shows the components it inherits from and which members it adds to the structs.

Because this mechanism is highly idiosyncratic, we're still exploring other ways to modularize CarboKitten.

## Model run
A model run is nothing but a loop running the `step!` function `n` times. The `run_model` method is overloaded to provide multiple levels of abstraction. In its most basic form, the `for` loop is abstracted into a `do` block, so we can run:

```julia
state = M.initial_state(input)
run_model(Model{M}, input, state) do step, frame
    ...
end
```

Note that the `state` variable is mutable, and we can use this to store/handle/inspect intermediate model states.

``` {.julia #run-model}
"""
    run_model(f, ::Type{Model{M}}, input::AbstractInput) where M
    run_model(f, ::Type{Model{M}}, input::AbstractInput, state::AbstractState) where M

Run a model and send `Frame`s to callback `f`. The type parameter `M` should be a model,
being a module with both `initial_state!` and `step!` defined.

The second version with explicit state is used by `H5Writer` so we can perform an additional
action between creating the initial state and starting the model run (saving metadata to the
HDF5 file).
"""
run_model(f, ::Type{Model{M}}, input::AbstractInput) where M =
    run_model(f, Model{M}, input, M.initial_state(input))

function run_model(f, ::Type{Model{M}}, input::AbstractInput, state::AbstractState) where M
    logger = get_logger(input)
    with_logger(logger) do
        step! = M.step!(input)

        @progress for w = 1:n_steps(input)
            f(w, step!(state))
        end
    end

    return state
end
```

``` {.julia file=src/RunModel.jl}
module RunModel

import ..CarboKitten: n_steps, run_model, get_logger, Model, AbstractInput, AbstractState
using ProgressLogging
using Logging: with_logger

<<run-model>>

end
```
