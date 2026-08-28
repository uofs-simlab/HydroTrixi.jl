@muladd begin
#! format: noindent

@doc raw"""
    IndicatorTotalVariation(semi; variable)

Compute the element-local total variation of the dimensionless sensor selected by
`variable` for a one-dimensional LGL-DGSEM discretization. On each element ``k``, the
indicator is
```math
\eta_k(t) \coloneqq \sum_{i=0}^N \omega_i
\left|\sum_{j=0}^N D_{ij}v_{k,j}(t)\right|
\approx \int_{-1}^1
\left|\frac{\mathrm{d}}{\mathrm{d}\xi}v_k^N(\xi,t)\right|\mathrm{d}\xi,
```
where ``N`` is the polynomial degree, ``\omega_i`` are the LGL quadrature weights, and
``D_{ij}`` are entries of the reference differentiation matrix. Any desired normalization
must be included in `variable` itself.
"""
struct IndicatorTotalVariation{Variable, Cache} <: Trixi.AbstractIndicator
    variable::Variable
    cache::Cache
end

function IndicatorTotalVariation(semi::Trixi.AbstractSemidiscretization;
                                 variable)
    cache = Trixi.create_cache(IndicatorTotalVariation, semi)
    return IndicatorTotalVariation(variable, cache)
end

function Trixi.create_cache(::Type{IndicatorTotalVariation},
                            equations::Trixi.AbstractEquations{1},
                            basis::Trixi.LobattoLegendreBasis)
    RealT = real(basis)
    alpha = Vector{RealT}()
    nodal_values_threaded = [Vector{RealT}(undef, Trixi.nnodes(basis))
                             for _ in 1:Threads.maxthreadid()]
    return (; alpha, nodal_values_threaded)
end

function (indicator::IndicatorTotalVariation)(u::AbstractArray{<:Any, 3},
                                              mesh::Trixi.TreeMesh{1}, equations,
                                              dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                              cache; kwargs...)
    (; alpha, nodal_values_threaded) = indicator.cache
    (; derivative_matrix, weights) = dg.basis
    resize!(alpha, Trixi.nelements(dg, cache))

    variable = indicator.variable
    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        nodal_values = nodal_values_threaded[Threads.threadid()]

        for i in Trixi.eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, element)
            nodal_values[i] = variable(u_local, equations)
        end

        total_variation = zero(eltype(alpha))
        for i in Trixi.eachnode(dg)
            derivative = zero(eltype(alpha))
            for j in Trixi.eachnode(dg)
                derivative += derivative_matrix[i, j] * nodal_values[j]
            end
            total_variation += weights[i] * abs(derivative)
        end
        alpha[element] = total_variation
    end

    return alpha
end
end # @muladd
