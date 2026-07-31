# Let Trixi.jl wrap contiguous implicit state views without copying
@inline function Trixi.storage_type(::Type{<:SubArray{<:Any, 1, <:Array,
                                                      <:Tuple{<:AbstractUnitRange}, true}})
    return Array
end

# Remake an implicit semidiscretization with caches matching the destination element type
# in order to evaluate a residual with dual numbers for automatic differentiation
function remake_implicit_semidiscretization(dual_cache_input)
    du_ode, semi = dual_cache_input
    return Trixi.remake(semi; uEltype = eltype(du_ode))
end

# Cache one remade semidiscretization for each dual-number type needed to evaluate the
# residual Jacobian with ForwardDiff.jl
mutable struct RHSImplicitCache{Cache}
    dual_semidiscretizations::Cache

    function RHSImplicitCache()
        dual_semidiscretizations = GeneralLazyBufferCache(remake_implicit_semidiscretization)
        return new{typeof(dual_semidiscretizations)}(dual_semidiscretizations)
    end
end

# Reconstruct a SemidiscretizationImplicit with the requested spatial cache element type
function Trixi.remake(semi::SemidiscretizationImplicit;
                      uEltype = eltype(semi.semi_base.cache.elements))
    semi_base = Trixi.remake(semi.semi_base; uEltype = uEltype)
    return SemidiscretizationImplicit(semi_base, semi.operator_temporal,
                                      semi.passive_variables)
end

# Evaluate a residual with caches matching the destination element type
function (rhs::RHSImplicitCache)(du_ode, u_ode, semi::SemidiscretizationImplicit, t)
    destination_eltype = eltype(du_ode)
    cache_eltype = eltype(semi.semi_base.cache.elements)

    if destination_eltype === cache_eltype
        rhs_implicit!(du_ode, u_ode, semi, t)
        return nothing
    end

    dual_semi = rhs.dual_semidiscretizations[(du_ode, semi)]
    rhs_implicit!(du_ode, u_ode, dual_semi, t)
    return nothing
end

# Accept ODE right-hand sides without an implicit element-type cache
invalidate_rhs_implicit_cache!(rhs) = nothing

# Discard cached dual-number semidiscretizations after a mesh topology change
function invalidate_rhs_implicit_cache!(rhs::RHSImplicitCache)
    rhs.dual_semidiscretizations = GeneralLazyBufferCache(remake_implicit_semidiscretization)
    return nothing
end

# Find the cached callable through nested ODE function wrappers
function invalidate_rhs_implicit_cache!(ode_function::SciMLBase.ODEFunction)
    return invalidate_rhs_implicit_cache!(ode_function.f)
end
