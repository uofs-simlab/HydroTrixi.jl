@doc raw"""
    HydrologicProblemRichardsClosedColumn(; soil_model, domain, tspan, base_head, amplitude)

Return a one-dimensional Richards-equation redistribution problem on a closed column with
homogeneous Neumann boundary conditions at both ends. The prescribed numerical boundary
fluxes vanish, so the semi-discrete water mass is constant. With [`MixedForm`](@ref), this
linear invariant is also retained by the time integration method up to roundoff.

The problem provides the constitutive map ``\theta = \vartheta(\psi)`` required by the
mixed formulation. The pressure head is initialized as
```math
\psi(z, 0) = \psi_{\text{base}} + z + A \left( 1 - \cos\left(\frac{2\pi (z - z_{\min})}{L}\right) \right),
```
where ``\psi_{\text{base}}`` is `base_head`, ``A`` is `amplitude`, and 
``L = z_{\max} - z_{\min}``. This profile satisfies ``\partial \psi / \partial z = 1`` at 
both boundaries, and therefore the flux
```math
f(\psi, \partial_z\psi) =
\mathcal{K}(\psi) \left( \frac{\partial \psi}{\partial z} - 1 \right)
```
vanishes at ``z = z_{\min}`` and ``z = z_{\max}`` at the initial time, and should
remain zero for all time due to the homogeneous Neumann boundary conditions.

The returned problem setup contains the fields `equations`, `state_to_evolved`,
`evolved_to_state`, `initial_condition`, `boundary_conditions`, `domain`, and `tspan`.
"""
function HydrologicProblemRichardsClosedColumn(; soil_model = default_soil_model(),
                                               domain = ((0.0,), (0.4,)),
                                               tspan = (0.0, 360.0), base_head = -1.0,
                                               amplitude = 0.05)
    lower, upper = domain
    z_min = lower[1]
    z_max = upper[1]
    length_scale = z_max - z_min
    if length_scale <= 0
        throw(ArgumentError("Expected a positive domain length."))
    end

    equations = RichardsEquation1D(soil_model = soil_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content

    function initial_condition(x, t, equations)
        z = x[1]
        phase = 2 * pi * (z - z_min) / length_scale
        psi = base_head + z + amplitude * (1 - cos(phase))
        return Trixi.SVector(psi)
    end

    zero_boundary_flux(x, t, equations) = Trixi.SVector(0.0)
    boundary_conditions = (; x_neg = Trixi.BoundaryConditionNeumann(zero_boundary_flux),
                           x_pos = Trixi.BoundaryConditionNeumann(zero_boundary_flux),)

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = initial_condition,
                             boundary_conditions = boundary_conditions, domain = domain,
                             tspan = tspan)
end
