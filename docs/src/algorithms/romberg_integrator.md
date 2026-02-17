Romberg Integration
===================

We use the `QuadGK.jl` package for quadrature. A fallback could be the following implementation of the Romberg method, which offers good compromise between performance and implementation complexity.

Implementation
--------------

From Wikipedia:

Using $h_n = \frac{(b-a)}{2^{n+1}}$, the method can be inductively defined by

$$\begin{align}
R(0,0) &= h_0 (f(a) + f(b)) \\
R(n,0) &= \tfrac{1}{2} R(n{-}1,\,0) + 2h_n \sum_{k=1}^{2^{n-1}} f(a + (2k-1)h_{n-1}) \\
R(n,m) &= R(n,\,m{-}1) + \tfrac{1}{4^m-1} (R(n,\,m{-}1) - R(n{-}1,\,m{-}1)) \\
&= \frac{1}{4^m-1} ( 4^m R(n,\,m{-}1) - R(n{-}1,\, m{-}1))
\end{align}$$

where $n \ge m$ and $m \ge 1$. In big O notation, the error for $R(n, m)$ is $O{\left(h_n^{2m+2}\right)}$.

``` {.julia file=src/Algorithms/Romberg.jl}
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
```
