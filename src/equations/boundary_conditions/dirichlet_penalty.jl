@muladd begin
#! format: noindent

@doc raw"""
    BoundaryConditionDirichletPenalty(boundary_value_function; penalty=1.0)

Dirichlet boundary condition with a configurable boundary penalty of the form
``\tau (u_{\text{inner}} - u_b)`` added to the boundary flux on the divergence pass for
parabolic diffusion operators. The `boundary_value_function` is called as
`boundary_value_function(x, t, equations)`, and `penalty` may be either a scalar or
a function called as `penalty(x, t, equations)`.

# References
- Arnold, D. N. (1982). An interior penalty finite element method with
  discontinuous elements. *SIAM Journal on Numerical Analysis*, 19(4), 742-760.
  [DOI: 10.1137/0719052](https://doi.org/10.1137/0719052)
- Arnold, D. N., Brezzi, F., Cockburn, B., Marini, L. D. (2002). Unified
  analysis of discontinuous Galerkin methods for elliptic problems.
  *SIAM Journal on Numerical Analysis*, 39(5), 1749-1779.
  [DOI: 10.1137/S0036142901384162](https://doi.org/10.1137/S0036142901384162)
"""
struct BoundaryConditionDirichletPenalty{BoundaryValueFunction, PenaltyFunction}
    boundary_value_function::BoundaryValueFunction
    penalty_function::PenaltyFunction
end

function BoundaryConditionDirichletPenalty(boundary_value_function; penalty = 1.0)
    penalty_function = penalty isa Function ? penalty : (x, t, equations) -> penalty
    return BoundaryConditionDirichletPenalty(boundary_value_function, penalty_function)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         orientation,
                                                                         direction,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         orientation,
                                                                         direction,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    penalty_strength = boundary_condition.penalty_function(x, t, equations)
    normal_sign = iseven(direction) ? one(penalty_strength) : -one(penalty_strength)
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         normal::AbstractVector,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         normal::AbstractVector,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    penalty_strength = boundary_condition.penalty_function(x, t, equations)
    normal_sign = sign(sum(normal))
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end
end # @muladd
