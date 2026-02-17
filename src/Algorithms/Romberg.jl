# ~/~ begin <<docs/src/algorithms/romberg_integrator.md#src/Algorithms/Romberg.jl>>[init]
module Romberg

"""
    romberg(::Type{ReturnType}, max_steps)(f, a, b, acc)
    
Implementation of the Romberg integrator, following a Python implementation on
Wikipedia. Reserves memory for refinement steps and returns a closure.

Function `f` is integrated over the domain `a` to `b` for a maximum amount of
refinement steps, or until accuracy `acc` is reached.
"""
function romberg(::Type{RA}, max_steps) where {RA}
    R1 = zeros(RA, max_steps)
    R2 = zeros(RA, max_steps)

    function (f, a::RB, b::RB, acc::RA) where {RB}
        RY = Base.return_types(f, (RB,))[1]
        
        Rp, Rc = R1, R2
        h::RB = (b - a) / 2
        ep::Int = 1
        Rp[1] = h * (f(a) + f(b))

        for i in 2:max_steps
            c::RY = zero(RY)

            for j in 1:ep
                c += f(a + (2 * j - 1) * h)
            end
            
            Rc[1] = h * c + 0.5 * Rp[1]

            n_k::Int = 1
            for j in 2:i
                n_k *= 4
                Rc[j] = (n_k * Rc[j - 1] - Rp[j - 1]) / (n_k - 1)
            end

            if i > 2 && abs(Rp[i-1] - Rc[i]) < acc
                return Rc[i]
            end

            ep *= 2
            h /= 2
            Rp, Rc = Rc, Rp
        end
        return Rp[max_steps]
    end
end

end  # module Romberg
# ~/~ end
