@muladd begin
#! format: noindent

@doc raw"""
    RichardsEquation1D(; soil_model)

A one-dimensional Richards equation model written in terms of the pressure head ``\psi``
for use with Trixi.jl's parabolic spatial discretization infrastructure. It represents the
spatial operator in either the pressure head form
```math
\frac{\partial \theta}{\partial \psi} \frac{\partial \psi}{\partial t} =
\frac{\partial}{\partial z}
\left( K(\psi) \left( \frac{\partial \psi}{\partial z} - 1 \right) \right),
```
or the mixed form
```math
\frac{\partial \theta}{\partial t} =
\frac{\partial}{\partial z}
\left( K(\psi) \left( \frac{\partial \psi}{\partial z} - 1 \right) \right),
```
where the constitutive relation between ``\theta`` and ``\psi`` is supplied by a separate
semidiscretization wrapper. The hydraulic conductivity is supplied through
`soil_model`, with `hydraulic_conductivity(psi, equations)` dispatching
on the model type parameter `SoilModel`. If `soil_model` is omitted, it
defaults to a [`Haverkamp`](@ref) model parameterized with the Celia (1990)
reference values reported in Ireson et al. (2023), Eq. (25), in SI units
(lengths in m, time in seconds).
"""
struct RichardsEquation1D{SoilModel} <:
       Trixi.AbstractEquationsParabolic{1, 1, Trixi.GradientVariablesConservative}
    soil_model::SoilModel
end

@inline default_soil_model() = Haverkamp(saturated_hydraulic_conductivity = 9.44e-5,
                                         alpha = 0.01936848004,
                                         beta = 3.96,
                                         A = 3.890790677e-4,
                                         gamma = 4.74,
                                         theta_s = 0.287,
                                         theta_r = 0.075)

@inline pressure_head(psi::Number) = psi
@inline pressure_head(u) = u[1]

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

@inline function Trixi.max_diffusivity(u, equations::RichardsEquation1D)
    return hydraulic_conductivity(u, equations.soil_model)
end

@doc raw"""
    water_content(u, equations::RichardsEquation1D)

Return the volumetric water content associated with the pressure head state `u` under the
Richards equation model `equations`.
"""
@inline function water_content(u, equations::RichardsEquation1D)
    soil_model = equations.soil_model
    S_e = effective_saturation(pressure_head(u), soil_model)
    return soil_model.theta_r + (soil_model.theta_s - soil_model.theta_r) * S_e
end

@inline function Trixi.flux(u, gradients, orientation::Integer,
                            equations::RichardsEquation1D)
    psi = u[1]
    dpsi_dz = first(gradients)[1]
    K_s = hydraulic_conductivity(psi, equations.soil_model)
    return K_s * (dpsi_dz - one(dpsi_dz))
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
