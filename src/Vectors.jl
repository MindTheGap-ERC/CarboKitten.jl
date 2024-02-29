# ~/~ begin <<docs/src/transport.md#src/Vectors.jl>>[init]
module Vectors

export Vec2

Vec2 = @NamedTuple{x::Float64, y::Float64}
Base.:+(a::Vec2, b::Vec2) = (x=a.x+b.x, y=a.y+b.y)
Base.abs2(a::Vec2) = a.x^2 + a.y^2
Base.abs(a::Vec2) = √(abs2(a))
Base.:*(a::Vec2, b::Float64) = (x=a.x*b, y=a.y*b)
Base.:/(a::Vec2, b::Float64) = (x=a.x/b, y=a.y/b)
Base.:*(a::Float64, b::Vec2) = b*a
Base.:-(a::Vec2, b::Vec2) = (x=a.x-b.x, y=a.y-b.y)
Base.:-(a::Vec2) = (x=-a.x, y=-a.y)
Base.zero(::Type{Vec2}) = (x=0.0, y=0.0)

end
# ~/~ end