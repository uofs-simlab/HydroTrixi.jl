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

# Source term corresponding to the manufactured pressure head
@inline function source_terms_richards_manufactured_solution(u, gradients, x, t,
                                                             equations)
    soil_model = equations.soil_model

    psi, psi_t, psi_z, psi_zz = richards_manufactured_profile(x, t)
    conductivity = hydraulic_conductivity(psi, soil_model)

    # Differentiate the Haverkamp conductivity for the manufactured source term
    if psi >= zero(psi)
        conductivity_derivative = zero(psi)
    else
        abs_psi = abs(psi)
        conductivity_denominator = soil_model.A + abs_psi^soil_model.gamma
        conductivity_derivative = soil_model.saturated_hydraulic_conductivity *
                                  soil_model.A * soil_model.gamma *
                                  abs_psi^(soil_model.gamma - 1) /
                                  conductivity_denominator^2
    end
    flux_derivative = conductivity_derivative * psi_z * (psi_z - 1) +
                      conductivity * psi_zz
    storage_derivative = water_capacity(psi, equations) * psi_t
    return Trixi.SVector(storage_derivative - flux_derivative)
end

@doc raw"""
    HydrologicProblemRichardsManufacturedSolution(; tspan = (0.0, 120.0),
                                                    soil_model = default_soil_model(),
                                                    penalty_factor = 1)

Return a one-dimensional manufactured-solution problem for the Richards equation. The
manufactured pressure head is
```math
\psi(z, t) =
0.204 \tanh\left(\frac{1}{2}\left(100z + \frac{t}{12} - 15\right)\right) - 0.411.
```
The boundary conditions impose this profile using the penalty formulation. The storage
source term
```math
s = C(\psi)\frac{\partial \psi}{\partial t}
- \frac{\partial}{\partial z}
\left(\mathcal{K}(\psi)\left(\frac{\partial \psi}{\partial z} - 1\right)\right),
\qquad C(\psi) = \vartheta'(\psi),
```
makes the profile solve the one-dimensional Richards equation with the Haverkamp
constitutive laws. The default setup uses the same Haverkamp parameters as
[`HydrologicProblemCelia1990`](@ref). The dimensionless `penalty_factor` is the coefficient
``C_\tau`` in the boundary penalty; setting it to zero omits the additional divergence-flux
penalty.

The problem uses depth ``z`` in metres on ``z \in [0, 0.2]`` and time in seconds on
``t \in [0, 120]`` by default. It is intended for regression and convergence checks of
mixed and pressure-head forms of the Richards equation.

# References
- List, F., Radu, F. A. (2016). A study on iterative methods for solving the
  Richards equation. *Computational Geosciences*, 20, 341-353.
  [DOI: 10.1007/s10596-016-9566-3](https://doi.org/10.1007/s10596-016-9566-3)
- Keita, S., Beljadid, A., Bourgault, Y. (2021). Implicit and semi-implicit
  second-order time stepping methods for the Richards equation.
  [arXiv:2105.05224](https://arxiv.org/abs/2105.05224)
"""
function HydrologicProblemRichardsManufacturedSolution(; tspan = (0.0, 120.0),
                                                       soil_model = default_soil_model(),
                                                       penalty_factor = 1)
    equations = RichardsEquation1D(soil_model = soil_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content
    boundary_condition = BoundaryConditionDirichletPenalty(richards_manufactured_solution;
                                                           penalty_factor)
    boundary_conditions = (; x_neg = boundary_condition, x_pos = boundary_condition)

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = richards_manufactured_solution,
                             boundary_conditions = boundary_conditions,
                             source_terms = source_terms_richards_manufactured_solution,
                             domain = ((0.0,), (0.2,)), tspan = tspan)
end
end # @muladd
