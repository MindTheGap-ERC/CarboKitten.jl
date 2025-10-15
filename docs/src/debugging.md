# Trouble shooting

As with any numerical model, there are many ways to shoot yourself in the foot using Entangled. Changing input paramaters may lead to numerical instablities and eventually even models crashing. When this happens it can be hard to understand what parameter settings or combinations thereof are responsible. Some tactics that have helped us in the past:

- Did you set a lithification or cementation time? The default cementation time is set to infinity (meaning the entire contents of the active layer is deposited on every iteration) which can lead to problems. Try a setting that is close in scale to your time step.
- Try reducing your time step. If you still encounter problems, but later in the simulation, it means that something in the settings is systematically driving the platform to numerical instability.
- If you supplied an additional (wave induced) velocity component to the transport model, this can lead to sediment build-up near the shore. In this case increasing the diffusivity will help to stabilize your model.

## Diagnostic mode

```component-dag
CarboKitten.Components.Production
```

You can enable diagnostic mode by setting `diagnostics = true` in the input settings. This will enable all debug messages, but some of these diagnostics can be expensive to compute.

## Implementation

``` {.julia file=src/Components/Diagnostics.jl}
@compose module Diagnostics

import ...CarboKitten: get_logger
using ..Common
using Logging
using TerminalLoggers: TerminalLogger
using LoggingExtras: TeeLogger, MinLevelLogger, FileLogger

@kwdef struct Input <: AbstractInput
    diagnostics::Bool = false
    log_file::String = "./carbokitten.log"
end

function get_logger(input::AbstractInput)
    if input.diagnostics
        io = open(input.log_file, "w+")
        file_logger = MinLevelLogger(FileLogger(input.log_file), Logging.Debug)
        terminal_logger = TerminalLogger(right_justify=80)
        return TeeLogger(file_logger, terminal_logger)
    else
        return TerminalLogger(right_justify=80)
    end
end

end
```
