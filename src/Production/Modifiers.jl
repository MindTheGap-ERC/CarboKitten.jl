# ~/~ begin <<docs/src/production/modifiers.md#src/Production/Modifiers.jl>>[init]
module Modifiers
    import ..Abstract: AbstractProduction, production_profile

    # ~/~ begin <<docs/src/production/modifiers.md#multiply-production>>[init]
    # =============================================================================
    # Time-window modifier — AbstractProduction transformer
    # =============================================================================
    
    const _ProdTime     = typeof(1.0u"Myr")
    const _ProdTimeSpec = Union{Colon, Tuple{_ProdTime,_ProdTime}}
    
    """
        MultiplyProduction(base, factor; t_range=:)
    
    Wraps `base::AbstractProduction`, multiplying its output by `factor` during
    `t_range`. Outside `t_range` the base production is unchanged.
    
    This implements the modifier pattern as `AbstractProduction -> AbstractProduction`:
    modifiers compose directly in the production spec rather than in a separate
    `production_modifiers` list on `Input`.
    """
    @kwdef struct MultiplyProduction <: AbstractProduction
        base::AbstractProduction
        factor::Float64
        t_range::_ProdTimeSpec = (:)
    end
    
    MultiplyProduction(base, factor::Real; kwargs...) =
        MultiplyProduction(; base=base, factor=Float64(factor), kwargs...)
    
    is_benthic(p::MultiplyProduction)      = is_benthic(p.base)
    is_pelagic(p::MultiplyProduction)      = is_pelagic(p.base)
    is_interpolated(p::MultiplyProduction) = is_interpolated(p.base)
    
    # ~/~ begin <<docs/src/production/modifiers.md#multiply-production-profile>>[init]
    function production_profile(input::AbstractInput, p::MultiplyProduction)
        base_profile = production_profile(input, p.base)
        return function(t, w)
            f = p.t_range isa Colon || (p.t_range[1] <= t <= p.t_range[2]) ? p.factor : 1.0
            return base_profile(t, w) * f
        end
    end
    # ~/~ end
    # ~/~ end
    # ~/~ begin <<docs/src/production/modifiers.md#multiply-production-profile>>[init]
    function production_profile(input::AbstractInput, p::MultiplyProduction)
        base_profile = production_profile(input, p.base)
        return function(t, w)
            f = p.t_range isa Colon || (p.t_range[1] <= t <= p.t_range[2]) ? p.factor : 1.0
            return base_profile(t, w) * f
        end
    end
    # ~/~ end
    # ~/~ begin <<docs/src/production/modifiers.md#production-boost>>[init]
    @kwdef struct ProductionBoost
        factor::Float64
        t_range::_ProdTimeSpec = (:)
    end
    
    Base.:*(p::AbstractProduction, b::ProductionBoost) = MultiplyProduction(p, b.factor, b.t_range)
    # ~/~ end
end
# ~/~ end
