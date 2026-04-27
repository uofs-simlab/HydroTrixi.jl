@muladd begin
#! format: noindent

@doc raw"""
    RichardsEquation1D(; hydraulic_conductivity_model)

A one-dimensional Richards equation model written in terms of the pressure head ``\psi``
for use with Trixi.jl's parabolic spatial discretization infrastructure. It represents the 
spatial operator in either the pressure head form 
```math
\frac{\partial \theta}{\partial \psi}\partial_t \psi = \partial_z \left( K(\psi) (\partial_z \psi - 1) \right),
```
or the mixed form
```math
\partial_t \theta = \partial_z \left( K(\psi) (\partial_z \psi - 1) \right),
```
where the constitutive relation between ``\theta`` and ``\psi`` is supplied by a separate 
semidiscretization wrapper. The hydraulic conductivity is supplied through
`hydraulic_conductivity_model`, with `hydraulic_conductivity(psi, equations)` dispatching
on the model type parameter `HydraulicConductivityModel`.
"""
struct RichardsEquation1D{HydraulicConductivityModel} <:
       Trixi.AbstractEquationsParabolic{1, 1, Trixi.GradientVariablesConservative}
    hydraulic_conductivity_model::HydraulicConductivityModel
end

@doc raw"""
    VanGenuchten(; saturated_hydraulic_conductivity, alpha, n,
                   m = 1 - 1 / n, pore_connectivity = 1 / 2)

A van Genuchten-Mualem hydraulic conductivity model. The effective saturation is
```math
S_e(\psi) =
\begin{cases}
1, & \psi \ge 0, \\
\left(1 + (\alpha |\psi|)^n\right)^{-m}, & \psi < 0,
\end{cases}
```
and the hydraulic conductivity is
```math
K(\psi) = K_s S_e(\psi)^l
\left(1 - \left(1 - S_e(\psi)^{1 / m}\right)^m\right)^2,
```
where ``K_s`` is the `saturated_hydraulic_conductivity` parameter and ``l`` is
the `pore_connectivity` parameter.
"""
struct VanGenuchten{RealT}
    saturated_hydraulic_conductivity::RealT
    alpha::RealT
    n::RealT
    m::RealT
    pore_connectivity::RealT
end

function VanGenuchten(; saturated_hydraulic_conductivity,
                      alpha,
                      n,
                      m = 1 - 1 / n,
                      pore_connectivity = 1 / 2)
    RealT = promote_type(typeof(saturated_hydraulic_conductivity), typeof(alpha),
                         typeof(n), typeof(m), typeof(pore_connectivity))
    return VanGenuchten{RealT}(saturated_hydraulic_conductivity, alpha, n, m,
                               pore_connectivity)
end

@inline _default_hydraulic_conductivity_model() = VanGenuchten(saturated_hydraulic_conductivity = 1.0,
                                                               alpha = 1.0,
                                                               n = 2.0)

function RichardsEquation1D(;
                            hydraulic_conductivity_model = _default_hydraulic_conductivity_model())
    return RichardsEquation1D(hydraulic_conductivity_model)
end

@inline Trixi.varnames(::typeof(Trixi.cons2cons), ::RichardsEquation1D) = ("psi",)
@inline Trixi.varnames(::typeof(Trixi.cons2prim), ::RichardsEquation1D) = ("psi",)
@inline Trixi.varnames(::typeof(Trixi.cons2entropy), ::RichardsEquation1D) = ("psi",)

@inline Trixi.cons2prim(u, ::RichardsEquation1D) = u
@inline Trixi.cons2entropy(u, ::RichardsEquation1D) = u

@inline function hydraulic_conductivity(psi,
                                        equations::RichardsEquation1D{HydraulicConductivityModel}) where {HydraulicConductivityModel}
    return hydraulic_conductivity(psi[1], equations.hydraulic_conductivity_model)
end

@inline function hydraulic_conductivity(psi, model::VanGenuchten)
    effective_saturation = if psi >= zero(psi)
        one(psi)
    else
        (one(psi) + (model.alpha * abs(psi))^model.n)^(-model.m)
    end
    conductivity_factor = (one(effective_saturation) -
                           (one(effective_saturation) -
                            effective_saturation^(inv(model.m)))^model.m)^2
    return model.saturated_hydraulic_conductivity *
           effective_saturation^model.pore_connectivity * conductivity_factor
end

@inline Trixi.have_constant_diffusivity(::RichardsEquation1D) = Trixi.False()

@inline function Trixi.max_diffusivity(u, equations::RichardsEquation1D)
    return hydraulic_conductivity(u, equations)
end

@inline function Trixi.flux(u, gradients, orientation::Integer,
                            equations::RichardsEquation1D)
    psi = u[1]
    dpsi_dz = first(gradients)[1]
    return hydraulic_conductivity(psi, equations) * (dpsi_dz - one(dpsi_dz))
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        normal::AbstractVector,
                                                                        x,
                                                                        t,
                                                                        operator_type::Trixi.Gradient,
                                                                        equations::RichardsEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        normal::AbstractVector,
                                                                        x,
                                                                        t,
                                                                        operator_type::Trixi.Divergence,
                                                                        equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        orientation,
                                                                        direction,
                                                                        x,
                                                                        t,
                                                                        operator_type::Trixi.Gradient,
                                                                        equations::RichardsEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        orientation,
                                                                        direction,
                                                                        x,
                                                                        t,
                                                                        operator_type::Trixi.Divergence,
                                                                        equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      normal::AbstractVector,
                                                                      x,
                                                                      t,
                                                                      operator_type::Trixi.Divergence,
                                                                      equations::RichardsEquation1D)
    return boundary_condition.boundary_normal_flux_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      normal::AbstractVector,
                                                                      x,
                                                                      t,
                                                                      operator_type::Trixi.Gradient,
                                                                      equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      orientation,
                                                                      direction,
                                                                      x,
                                                                      t,
                                                                      operator_type::Trixi.Divergence,
                                                                      equations::RichardsEquation1D)
    return boundary_condition.boundary_normal_flux_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      orientation,
                                                                      direction,
                                                                      x,
                                                                      t,
                                                                      operator_type::Trixi.Gradient,
                                                                      equations::RichardsEquation1D)
    return flux_inner
end
end # @muladd
