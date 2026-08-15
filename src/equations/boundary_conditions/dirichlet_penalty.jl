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
f^* = f_{\mathrm{inner}} - n_b C_\tau\,\kappa(u_b)\frac{N(N+1)}{h}(u_{\mathrm{inner}} - u_b)
```
on the divergence pass, where ``n_b \in \{-1, 1\}`` is
the outward unit normal. The solver supplies the polynomial degree ``N`` and boundary-cell
size ``h``. Since the endpoint weight of the collocated Legendre-Gauss-Lobatto rule on
``[-1, 1]`` is ``2/(N(N+1))`` and the inverse Jacobian of the affine map from the reference 
cell to the physical boundary cell is ``2/h``, jumps in the numerical trace get lifted with 
a scaling of ``N(N+1)/h`` when converted into a gradient on the reference element, and 
further scaled by the diffusion coefficient to obtain a flux on the divergence pass, where 
``\kappa(u_b)`` is the diffusion coefficient evaluated at the boundary value using 
`boundary_penalty_coefficient(u_b, equations)`. The dimensionless `penalty_factor` is
``C_\tau``; its default value is one. Setting it to zero omits this additional
divergence-flux penalty while retaining the prescribed trace on the gradient pass.

# References
- Manzanero, J., Rueda-Ramírez, A. M., Rubio, G., Ferrer, E. (2018). The Bassi Rebay 1
  scheme is a special case of the symmetric interior penalty formulation for discontinuous
  Galerkin discretisations with Gauss-Lobatto points. *Journal of Computational Physics*,
  363, 1-10.
  [DOI: 10.1016/j.jcp.2018.02.035](https://doi.org/10.1016/j.jcp.2018.02.035)
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
