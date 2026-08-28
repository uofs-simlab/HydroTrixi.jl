@muladd begin
#! format: noindent

@doc raw"""
    RichardsEquation1D(; soil_model)

A one-dimensional Richards equation model for vertical flow in a soil column, with depth
``z`` measured positive downward. The pressure head is ``\psi(z,t)``, the constitutive
water-content function is ``\vartheta(\psi)``, and the nonlinear flux is
``f(\psi, \partial_z\psi) \coloneqq \kappa(\psi)(\partial_z\psi - 1)``. The unit
gravitational-gradient term follows from the downward-positive depth convention.

The model supplies the spatial operator shared by the pressure-head formulation
```math
c(\psi) \partial_t \psi =
\partial_z \left( \kappa(\psi) \left( \partial_z \psi - 1 \right) \right),
```
where ``c(\psi) \coloneqq \vartheta'(\psi)``, and the mixed formulation
```math
\partial_t \theta =
\partial_z \left( \kappa(\psi) \left( \partial_z \psi - 1 \right) \right),
\qquad
\theta = \vartheta(\psi).
```
The temporal formulation and constitutive constraint are supplied by
[`SemidiscretizationImplicit`](@ref). The hydraulic conductivity is supplied through
`soil_model`, with `hydraulic_conductivity(psi, equations)` dispatching on the model type
parameter `SoilModel`. [`BoundaryConditionDirichletPenalty`](@ref) uses the default penalty
``\kappa(\psi_{\mathrm{D}})N(N+1)/h``, where ``\psi_{\mathrm{D}}`` is the prescribed
boundary pressure head. If `soil_model` is omitted, it defaults to a [`Haverkamp`](@ref)
model parameterized with the Celia et al. (1990) reference values reported in Ireson et al.
(2023), Eq. (25), in SI units (lengths in metres and time in seconds).
"""
struct RichardsEquation1D{SoilModel} <:
       Trixi.AbstractEquationsParabolic{1, 1, Trixi.GradientVariablesConservative}
    soil_model::SoilModel
end

@inline default_soil_model() = Haverkamp(saturated_hydraulic_conductivity = 9.44e-5,
                                         alpha = 0.01936848004, beta = 3.96,
                                         A = 3.890790677e-4, gamma = 4.74,
                                         theta_s = 0.287,
                                         theta_r = 0.075)

@doc raw"""
    pressure_head(u)
    pressure_head(u, equations::RichardsEquation1D)

Return the pressure head ``\psi`` stored in state `u`.

For scalar states, `u` is returned directly. For vector-like states, the first component is
interpreted as pressure head. For example, the mixed formulation stores pressure head in the
first component and water content in the second component, so `pressure_head(u)` returns the
first component.
"""
@inline pressure_head(psi::Number) = psi
@inline pressure_head(u) = u[1]
@inline pressure_head(u, ::RichardsEquation1D) = pressure_head(u)

@inline function pressure_head_out_of_domain(u, semi, t, ::RichardsEquation1D)
    pressure_heads = state_variable_view(u, semi)
    return any(psi -> psi >= zero(psi), pressure_heads)
end

function RichardsEquation1D(; soil_model = default_soil_model())
    return RichardsEquation1D(soil_model)
end

@inline Trixi.varnames(::typeof(Trixi.cons2cons), ::RichardsEquation1D) = ("psi",)
@inline Trixi.varnames(::typeof(Trixi.cons2prim), ::RichardsEquation1D) = ("psi",)
@inline Trixi.varnames(::typeof(Trixi.cons2entropy), ::RichardsEquation1D) = ("psi",)
@inline Trixi.default_analysis_integrals(::RichardsEquation1D) = ()

@inline Trixi.cons2prim(u, ::RichardsEquation1D) = u
@inline Trixi.cons2entropy(u, ::RichardsEquation1D) = u
@inline Trixi.have_constant_diffusivity(::RichardsEquation1D) = Trixi.False()

@doc raw"""
    hydraulic_conductivity(u, equations::RichardsEquation1D)

Return the hydraulic conductivity associated with the pressure head stored in `u`.
"""
@inline function hydraulic_conductivity(u, equations::RichardsEquation1D)
    return hydraulic_conductivity(pressure_head(u), equations.soil_model)
end

@inline function Trixi.max_diffusivity(u, equations::RichardsEquation1D)
    return hydraulic_conductivity(u, equations)
end

@doc raw"""
    effective_saturation(u, equations::RichardsEquation1D)

Return the effective saturation associated with the pressure head ``\psi`` stored in `u`,
defined as
```math
S_{\mathrm{e}}(\psi) \coloneqq
\frac{\vartheta(\psi) - \theta_{\mathrm{r}}}
     {\theta_{\mathrm{s}} - \theta_{\mathrm{r}}},
```
where ``\theta_{\mathrm{r}}`` and ``\theta_{\mathrm{s}}`` are the residual and saturated
water-content values, respectively, for the Richards equation model `equations`.
"""
@inline function effective_saturation(u, equations::RichardsEquation1D)
    return effective_saturation(pressure_head(u), equations.soil_model)
end

@doc raw"""
    water_content(u, equations::RichardsEquation1D)

Return the volumetric water content ``\vartheta(\psi)`` associated with the pressure head
state `u` under the Richards equation model `equations`.
"""
@inline function water_content(u, equations::RichardsEquation1D)
    soil_model = equations.soil_model
    S_e = effective_saturation(u, equations)
    return soil_model.theta_r + (soil_model.theta_s - soil_model.theta_r) * S_e
end

@doc raw"""
    water_capacity(u, equations::RichardsEquation1D)

Return the capacity ``c(\psi) \coloneqq \vartheta'(\psi)`` associated with the pressure
head state `u` under the Richards equation model `equations`.
"""
@inline function water_capacity(u, equations::RichardsEquation1D)
    return water_capacity(pressure_head(u), equations.soil_model)
end

@inline function water_capacity(psi, model::Haverkamp)
    if psi >= zero(psi)
        return zero(psi)
    end

    abs_psi = abs(psi)
    theta_range = model.theta_s - model.theta_r
    denominator = model.alpha + abs_psi^model.beta
    return theta_range * model.alpha * model.beta * abs_psi^(model.beta - 1) /
           denominator^2
end

@inline function water_capacity(psi, model::VanGenuchten)
    if psi >= zero(psi)
        return zero(psi)
    end

    abs_psi = abs(psi)
    theta_range = model.theta_s - model.theta_r
    saturation_denominator = one(psi) + (model.alpha * abs_psi)^model.n
    return theta_range * model.m * model.n * model.alpha^model.n *
           abs_psi^(model.n - 1) * saturation_denominator^(-model.m - 1)
end

@doc raw"""
    pressure_head_from_water_content(theta, equations::RichardsEquation1D)

Return the pressure head associated with water content `theta` for the retention curve
stored in `equations`.
"""
@inline function pressure_head_from_water_content(theta, equations::RichardsEquation1D)
    return pressure_head_from_water_content(theta, equations.soil_model)
end

@inline function pressure_head_from_water_content(theta, model::Haverkamp)
    effective_saturation = (theta - model.theta_r) / (model.theta_s - model.theta_r)
    unsaturated_head = (model.alpha *
                        (one(effective_saturation) - effective_saturation) /
                        effective_saturation)^(inv(model.beta))
    return -unsaturated_head
end

@inline function Trixi.flux(u, gradients, orientation::Integer,
                            equations::RichardsEquation1D)
    psi = u[1]
    dpsi_dz = first(gradients)[1]
    K_s = hydraulic_conductivity(psi, equations.soil_model)
    return K_s * (dpsi_dz - one(dpsi_dz))
end

# Penalty coefficient used in BoundaryConditionDirichletPenalty
@inline function boundary_penalty_coefficient(u_boundary, equations::RichardsEquation1D)
    return hydraulic_conductivity(u_boundary, equations)
end

# These boundary condition methods are needed because Trixi.jl defines the corresponding
# generic behavior only for `AbstractLaplaceDiffusion`, which is not a supertype of
# `RichardsEquation1D`. The behavior here matches those Laplace-diffusion methods.
@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        normal::AbstractVector,
                                                                        x, t,
                                                                        operator_type::Trixi.Gradient,
                                                                        equations::RichardsEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        normal::AbstractVector,
                                                                        x, t,
                                                                        operator_type::Trixi.Divergence,
                                                                        equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        orientation,
                                                                        direction, x, t,
                                                                        operator_type::Trixi.Gradient,
                                                                        equations::RichardsEquation1D)
    return boundary_condition.boundary_value_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionDirichlet)(flux_inner,
                                                                        u_inner,
                                                                        orientation,
                                                                        direction, x, t,
                                                                        operator_type::Trixi.Divergence,
                                                                        equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      normal::AbstractVector,
                                                                      x, t,
                                                                      operator_type::Trixi.Divergence,
                                                                      equations::RichardsEquation1D)
    return boundary_condition.boundary_normal_flux_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      normal::AbstractVector,
                                                                      x, t,
                                                                      operator_type::Trixi.Gradient,
                                                                      equations::RichardsEquation1D)
    return flux_inner
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      orientation,
                                                                      direction, x, t,
                                                                      operator_type::Trixi.Divergence,
                                                                      equations::RichardsEquation1D)
    return boundary_condition.boundary_normal_flux_function(x, t, equations)
end

@inline function (boundary_condition::Trixi.BoundaryConditionNeumann)(flux_inner,
                                                                      u_inner,
                                                                      orientation,
                                                                      direction, x, t,
                                                                      operator_type::Trixi.Gradient,
                                                                      equations::RichardsEquation1D)
    return flux_inner
end
end # @muladd
