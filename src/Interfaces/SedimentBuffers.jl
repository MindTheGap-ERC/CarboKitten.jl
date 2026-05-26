# ~/~ begin <<docs/src/algorithms/stack_interface.md#src/Interfaces/SedimentBuffers.jl>>[init]
module SedimentBuffers

export push_sediment!, pop_sediment!, peek_sediment
using ..Chunks

"""
    push_sediment!(buffer, sediment)
    push_sediment!(buffer, sediment, chunk)

Push deposition onto the sediment buffer.
"""
function push_sediment! end

push_sediment!(buffer, sediment) = push_sediment!(buffer, sediment, Serial(Chunk(:, :)))

"""
    pop_sediment!(buffer, amount, out)
    pop_sediment!(buffer, amount, out, chunk)

Pop amount off the buffer, writing to sediment array.
"""
function pop_sediment! end

pop_sediment!(buffer, amount, out) = pop_sediment(buffer, amount, out, Serial(Chunk(:, :)))

"""
    peek_sediment(buffer, amount, out)
    peek_sediment(buffer, amount, out, chunk)

Peek at the top of the buffer.
"""
function peek_sediment end

peek_sediment(buffer, amount, out) = peek_sediment(buffer, amount, out, Serial(Chunk(:, :)))

end
# ~/~ end
