@muladd begin
#! format: noindent

"""
    BoundaryConditionDirichletPenalty(boundary_value_function; penalty=1.0)

Dirichlet boundary condition with a configurable boundary penalty for parabolic
diffusion operators. The `boundary_value_function` is called as
`boundary_value_function(x, t, equations)`.

`penalty` may be either a scalar or a function called as
`penalty(x, t, equations)`.
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
                                                                         equations::AbstractDiffusionEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         orientation,
                                                                         direction,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::AbstractDiffusionEquation1D)
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
                                                                         equations::AbstractDiffusionEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         normal::AbstractVector,
                                                                         x,
                                                                         t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::AbstractDiffusionEquation1D)
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    penalty_strength = boundary_condition.penalty_function(x, t, equations)
    normal_sign = sign(sum(normal))
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end
end # @muladd
