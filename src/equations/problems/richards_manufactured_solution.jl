@muladd begin
#! format: noindent

# Manufactured pressure head used by `HydrologicProblemRichardsManufacturedSolution`
function richards_manufactured_solution(x, t, equations)
    psi, _, _, _ = richards_manufactured_profile(x, t)
    return Trixi.SVector(psi)
end

# Manufactured profile and derivatives in pressure-head form
@inline function richards_manufactured_profile(x, t)
    z = x[1]
    head_amplitude = 0.204
    head_offset = -0.411
    eta_z = 50.0
    eta_t = 1 / 24
    eta = 0.5 * (100 * z + t / 12 - 15)
    tanh_eta = tanh(eta)
    sech2_eta = 1 - tanh_eta^2

    psi = head_amplitude * tanh_eta + head_offset
    psi_t = head_amplitude * sech2_eta * eta_t
    psi_z = head_amplitude * sech2_eta * eta_z
    psi_zz = -2 * head_amplitude * eta_z^2 * tanh_eta * sech2_eta
    return psi, psi_t, psi_z, psi_zz
end

# Exact normal flux at the right boundary for a manufactured Neumann condition
@inline function richards_manufactured_right_boundary_flux(x, t, equations)
    psi, _, psi_z, _ = richards_manufactured_profile(x, t)
    flux = hydraulic_conductivity(psi, equations) * (psi_z - one(psi_z))
    return Trixi.SVector(flux)
end

# Source term corresponding to the manufactured pressure head
@inline function source_terms_richards_manufactured_solution(u, gradients, x, t,
                                                             equations)
    constitutive_model = equations.constitutive_model

    psi, psi_t, psi_z, psi_zz = richards_manufactured_profile(x, t)
    conductivity = hydraulic_conductivity(psi, constitutive_model)

    # Differentiate the Haverkamp conductivity for the manufactured source term
    if psi >= zero(psi)
        conductivity_derivative = zero(psi)
    else
        abs_psi = abs(psi)
        scaled_abs_psi = constitutive_model.b * abs_psi
        conductivity_denominator = one(psi) + scaled_abs_psi^constitutive_model.gamma
        conductivity_derivative = constitutive_model.saturated_hydraulic_conductivity *
                                  constitutive_model.b * constitutive_model.gamma *
                                  scaled_abs_psi^(constitutive_model.gamma - 1) /
                                  conductivity_denominator^2
    end
    flux_derivative = conductivity_derivative * psi_z * (psi_z - 1) +
                      conductivity * psi_zz
    storage_derivative = water_capacity(psi, equations) * psi_t
    return Trixi.SVector(storage_derivative - flux_derivative)
end

@doc raw"""
    HydrologicProblemRichardsManufacturedSolution(; tspan = (0.0, 120.0),
                                                    constitutive_model = default_constitutive_model(),
                                                    penalty_factor = 1,
                                                    boundary_conditions = :dirichlet_dirichlet)

Return a one-dimensional manufactured-solution problem for the Richards equation. The
manufactured pressure head is
```math
\psi(z, t) =
0.204 \tanh\left(\frac{1}{2}\left(100z + \frac{t}{12} - 15\right)\right) - 0.411.
```
The default `boundary_conditions = :dirichlet_dirichlet` imposes this profile at both
boundaries using penalty Dirichlet conditions. With `:dirichlet_neumann`, the left
boundary uses the same penalty Dirichlet condition and the right boundary imposes the
exact normal flux. A custom `(; x_neg, x_pos)` tuple overrides both boundaries;
`nothing` selects the default pair.

```julia
problem = HydrologicProblemRichardsManufacturedSolution(
    boundary_conditions = :dirichlet_neumann)
exact_solution = problem.initial_condition
```
The manufactured solution is imposed through the source term
```math
s(z,t) = c(\psi(z,t))\partial_t \psi(z,t)
- \partial_z \left(\kappa(\psi(z,t))
\left(\partial_z \psi(z,t) - 1\right)\right),
\qquad c(\psi) \coloneqq \vartheta'(\psi),
```
The default setup uses the same Haverkamp parameters as
[`HydrologicProblemCeliaHaverkamp`](@ref). The dimensionless `penalty_factor` is the
coefficient ``C_\tau`` for the Dirichlet boundaries in either named choice; setting it to
zero omits the divergence flux penalty term. Custom boundary tuples retain their own
penalty settings.

The problem uses depth ``z`` in metres on ``z \in [0, 0.2]`` and time in seconds on
``t \in [0, 120]`` by default. It is intended for regression and convergence checks of
mixed and pressure-head forms of the Richards equation.

# References
- Keita, S., Beljadid, A., Bourgault, Y. (2021). Implicit and semi-implicit
  second-order time stepping methods for the Richards equation.
  [arXiv:2105.05224](https://arxiv.org/abs/2105.05224)
"""
function HydrologicProblemRichardsManufacturedSolution(; tspan = (0.0, 120.0),
                                                       constitutive_model = default_constitutive_model(),
                                                       penalty_factor = 1,
                                                       boundary_conditions = :dirichlet_dirichlet)
    equations = RichardsEquation1D(constitutive_model = constitutive_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content
    if isnothing(boundary_conditions) || boundary_conditions isa Symbol
        if isnothing(boundary_conditions)
            boundary_conditions = :dirichlet_dirichlet
        end
        left_boundary = BoundaryConditionDirichletPenalty(richards_manufactured_solution;
                                                          penalty_factor)
        right_boundary = if boundary_conditions === :dirichlet_dirichlet
            left_boundary
        elseif boundary_conditions === :dirichlet_neumann
            Trixi.BoundaryConditionNeumann(richards_manufactured_right_boundary_flux)
        else
            throw(ArgumentError("Unknown boundary choice $boundary_conditions; use " *
                                ":dirichlet_dirichlet or :dirichlet_neumann."))
        end
        boundary_conditions = (; x_neg = left_boundary, x_pos = right_boundary)
    end

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = richards_manufactured_solution,
                             boundary_conditions = boundary_conditions,
                             source_terms = source_terms_richards_manufactured_solution,
                             domain = ((0.0,), (0.2,)), tspan = tspan)
end
end # @muladd
