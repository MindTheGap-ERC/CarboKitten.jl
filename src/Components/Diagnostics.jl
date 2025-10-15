# ~/~ begin <<docs/src/debugging.md#src/Components/Diagnostics.jl>>[init]
@compose module Diagnostics

using ..Common
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
        terminal_logger = MinLevelLogger(TerminalLogger(right_justify=80), Logging.Info)
        return TeeLogger(file_logger, terminal_logger)
    else
        return TerminalLogger(right_justify=80)
    end
end

end
# ~/~ end
