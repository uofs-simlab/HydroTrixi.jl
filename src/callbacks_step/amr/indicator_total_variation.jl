@muladd begin
#! format: noindent

@doc raw"""
    IndicatorTotalVariation(semi; variable, normalize = false)

Compute the element-local total variation of the sensor selected by `variable`
for a one-dimensional LGL-DGSEM discretization. On each element ``k``, the
indicator is
```math
\eta_k(t) \coloneqq \sum_{i=0}^N \omega_i
\left|\sum_{j=0}^N D_{ij}v_{k,j}(t)\right|
\approx \int_{-1}^1
\left|\frac{\mathrm{d}}{\mathrm{d}\xi}v_k^N(\xi,t)\right|\mathrm{d}\xi,
```
where ``N`` is the polynomial degree, ``\omega_i`` are the LGL quadrature weights, and
``D_{ij}`` are entries of the reference differentiation matrix.

With `normalize = true`, divide by the range of `variable` over all mesh nodes,
recomputed on every indicator evaluation:
```math
\widehat\eta_k = \frac{\eta_k}{v_{\max}-v_{\min}+\operatorname{eps}}.
```
Here `eps` is the machine epsilon of the indicator's floating-point type.
Constant fields have zero indicator. Extrema use the current nodal solution over
all MPI ranks, without adding boundary values or interelement jumps.
The normalized indicator is dimensionless; with `normalize = false`, it has the
same units as `variable` and retains the unnormalized behavior.
"""
struct IndicatorTotalVariation{Variable, Cache} <: Trixi.AbstractIndicator
    variable ::Variable
    normalize::Bool
    cache    ::Cache
end

function IndicatorTotalVariation(semi::Trixi.AbstractSemidiscretization;
                                 variable, normalize::Bool = false)
    cache = Trixi.create_cache(IndicatorTotalVariation, semi)
    return IndicatorTotalVariation(variable, normalize, cache)
end

function Trixi.create_cache(::Type{IndicatorTotalVariation},
                            equations::Trixi.AbstractEquations{1},
                            basis::Trixi.LobattoLegendreBasis)
    RealT = real(basis)
    alpha = Vector{RealT}()
    minima = Vector{RealT}()
    maxima = Vector{RealT}()
    nodal_values_threaded = [Vector{RealT}(undef, Trixi.nnodes(basis))
                             for _ in 1:Threads.maxthreadid()]
    return (; alpha, minima, maxima, nodal_values_threaded)
end

function (indicator::IndicatorTotalVariation)(u::AbstractArray{<:Any, 3},
                                              mesh::Trixi.TreeMesh{1}, equations,
                                              dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                              cache; kwargs...)
    (; variable, normalize) = indicator
    (; alpha, minima, maxima, nodal_values_threaded) = indicator.cache
    (; derivative_matrix, weights) = dg.basis
    resize!(alpha, Trixi.nelements(dg, cache))
    if normalize
        resize!(minima, length(alpha))
        resize!(maxima, length(alpha))
    end

    Trixi.@threaded for element in Trixi.eachelement(dg, cache)
        nodal_values = nodal_values_threaded[Threads.threadid()]

        for i in Trixi.eachnode(dg)
            u_local = Trixi.get_node_vars(u, equations, dg, i, element)
            nodal_values[i] = variable(u_local, equations)
        end

        if normalize
            minima[element], maxima[element] = extrema(nodal_values)
            # Differentiate differences so constant fields give exactly zero.
            offset = first(nodal_values)
            for i in Trixi.eachnode(dg)
                nodal_values[i] -= offset
            end
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

    if normalize
        v_min = minimum(minima; init = typemax(eltype(minima)))
        v_max = maximum(maxima; init = typemin(eltype(maxima)))
        if Trixi.mpi_isparallel()
            v_min = Trixi.MPI.Allreduce(v_min, min, Trixi.mpi_comm())
            v_max = Trixi.MPI.Allreduce(v_max, max, Trixi.mpi_comm())
        end
        scale = v_max - v_min + eps(eltype(alpha))
        alpha ./= scale
    end

    return alpha
end
end # @muladd
