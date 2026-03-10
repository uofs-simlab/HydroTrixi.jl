struct LinearDiffusionEquation1D{RealT <: Real} <: Trixi.AbstractLaplaceDiffusion{1, 1}
    diffusivity::RealT
end

Trixi.varnames(::typeof(Trixi.cons2cons), ::LinearDiffusionEquation1D) = ("scalar",)

@inline function Trixi.flux(u, gradients, orientation::Integer,
                            equations::LinearDiffusionEquation1D)
    dudx, = gradients
    return equations.diffusivity * dudx
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner, u_inner,
                                                                         orientation,
                                                                         direction,
                                                                         x, t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::LinearDiffusionEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner, u_inner,
                                                                         orientation,
                                                                         direction,
                                                                         x, t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::LinearDiffusionEquation1D)
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    penalty_strength = boundary_condition.penalty_function(x, t, equations)
    normal_sign = iseven(direction) ? one(penalty_strength) : -one(penalty_strength)
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner, u_inner,
                                                                         normal::AbstractVector,
                                                                         x, t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::LinearDiffusionEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner, u_inner,
                                                                         normal::AbstractVector,
                                                                         x, t,
                                                                         operator_type::Trixi.Divergence,
                                                                         equations::LinearDiffusionEquation1D)
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    penalty_strength = boundary_condition.penalty_function(x, t, equations)
    normal_sign = sign(sum(normal))
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end
