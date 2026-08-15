@muladd begin
#! format: noindent

@doc raw"""
    BoundaryConditionDirichletPenalty(boundary_value_function; penalty_factor=1)

Dirichlet boundary condition for parabolic operators with the numerical trace
```math
u^* = u_b
```
on the gradient pass and the numerical flux
```math
f^* = f_{\mathrm{inner}} - n_b C_\tau\,\kappa(u_b)\frac{(N+1)^2}{h}(u_{\mathrm{inner}} - u_b)
```
on the divergence pass, where ``n_b \in \{-1, 1\}`` is
the outward unit normal. The solver supplies the polynomial degree ``N`` and boundary-cell
size ``h``. The (possibly nonlinear) diffusion coefficient ``\kappa(u_b)`` is obtained from
`boundary_penalty_coefficient(u_b, equations)`. This returns the constant diffusivity for
linear diffusion and the hydraulic conductivity ``\mathcal{K}(\psi_D)`` for the Richards
equation. The specified `boundary_value_function` is called as
`boundary_value_function(x, t, equations)`. The dimensionless `penalty_factor` is
``C_\tau``; its default value is one, while zero recovers the unpenalized divergence flux.

# References
- Arnold, D. N. (1982). An interior penalty finite element method with
  discontinuous elements. *SIAM Journal on Numerical Analysis*, 19(4), 742-760.
  [DOI: 10.1137/0719052](https://doi.org/10.1137/0719052)
- Arnold, D. N., Brezzi, F., Cockburn, B., Marini, L. D. (2002). Unified
  analysis of discontinuous Galerkin methods for elliptic problems.
  *SIAM Journal on Numerical Analysis*, 39(5), 1749-1779.
  [DOI: 10.1137/S0036142901384162](https://doi.org/10.1137/S0036142901384162)
"""
struct BoundaryConditionDirichletPenalty{BoundaryValueFunction, PenaltyFactor}
    boundary_value_function::BoundaryValueFunction
    penalty_factor::PenaltyFactor
end

function BoundaryConditionDirichletPenalty(boundary_value_function;
                                           penalty_factor::Real = 1)
    if !isfinite(penalty_factor) || penalty_factor < zero(penalty_factor)
        throw(ArgumentError("`penalty_factor` must be finite and non-negative."))
    end
    return BoundaryConditionDirichletPenalty(boundary_value_function, penalty_factor)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         orientation,
                                                                         direction, x,
                                                                         t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichletPenalty)(flux_inner,
                                                                         u_inner,
                                                                         normal::AbstractVector,
                                                                         x, t,
                                                                         operator_type::Trixi.Gradient,
                                                                         equations::Trixi.AbstractEquationsParabolic{1})
    return boundary_condition.boundary_value_function(x, t, equations)
end

# Supply the equation-dependent coefficient in the penalty strength.
@inline function boundary_penalty_coefficient(u_boundary,
                                              equations::Trixi.AbstractLaplaceDiffusion)
    return Trixi.max_diffusivity(equations)
end

@inline function boundary_flux_divergence(boundary_condition::BoundaryConditionDirichletPenalty,
                                          flux_inner, u_inner, orientation, direction,
                                          x, t,
                                          equations::Trixi.AbstractEquationsParabolic{1},
                                          penalty_scale)
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)
    diffusion_coefficient = boundary_penalty_coefficient(u_boundary, equations)
    penalty_strength = boundary_condition.penalty_factor * diffusion_coefficient *
                       penalty_scale
    normal_sign = iseven(direction) ? one(penalty_strength) : -one(penalty_strength)
    return flux_inner - normal_sign * penalty_strength * (u_inner - u_boundary)
end
end # @muladd
