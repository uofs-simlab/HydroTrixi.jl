@doc raw"""
    HydrologicProblemCelia1990(; tspan = (0.0, 360.0))

Return the one-dimensional Richards-equation infiltration problem introduced by Celia,
Bouloutas, and Zarba (1990) as a `HydrologicProblem` for HydroTrixi.jl semidiscretizations.

The problem uses the Haverkamp constitutive laws in the form reported by Ireson et al.
(2023), Eq. (25):
```math
\vartheta(\psi) = \vartheta_r + (\vartheta_s - \vartheta_r)
\frac{\alpha}{\alpha + |\psi|^\beta},
\qquad
\mathcal{K}(\psi) = \mathcal{K}_s\,\frac{A}{A + |\psi|^\gamma},
```
for ``\psi < 0``, with saturated values ``\vartheta = \vartheta_s`` and
``\mathcal{K} = \mathcal{K}_s`` for ``\psi \ge 0``. The setup uses depth ``z`` in
metres, positive downward on ``z \in [0, 0.4]``, and time in seconds on
``t \in [0, 360]``. The pressure head is initialized at ``-0.615`` m, with a Dirichlet
boundary value of ``-0.207`` m at the soil surface (`x_neg`) and ``-0.615`` m at the
bottom of the column (`x_pos`).

The returned problem setup contains the fields `equations`, `state_to_evolved`,
`evolved_to_state`, `initial_condition`, `boundary_conditions`, `domain`, and `tspan`.

# References
- Celia, M. A., Bouloutas, E. T., Zarba, R. L. (1990). A general
  mass-conservative numerical solution for the unsaturated flow equation.
  *Water Resources Research*, 26(7), 1483-1496.
  [DOI: 10.1029/WR026i007p01483](https://doi.org/10.1029/WR026i007p01483)
- Ireson, A. M., Spiteri, R. J., Clark, M. P., Mathias, S. A. (2023).
  A simple, efficient, mass-conservative approach to solving Richards'
  equation (openRE, v1.0). *Geoscientific Model Development*, 16, 659-677.
  [DOI: 10.5194/gmd-16-659-2023](https://doi.org/10.5194/gmd-16-659-2023)
"""
function HydrologicProblemCelia1990(; tspan = (0.0, 360.0))
    soil_model = Haverkamp(saturated_hydraulic_conductivity = 9.44e-5,
                           alpha = 0.01936848004, beta = 3.96, A = 3.890790677e-4,
                           gamma = 4.74, theta_s = 0.287, theta_r = 0.075)
    equations = RichardsEquation1D(soil_model = soil_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content
    initial_condition(x, t, equations) = Trixi.SVector(-0.615)
    top_boundary_value(x, t, equations) = Trixi.SVector(-0.207)
    bottom_boundary_value(x, t, equations) = Trixi.SVector(-0.615)
    boundary_conditions = (; x_neg = Trixi.BoundaryConditionDirichlet(top_boundary_value),
                           x_pos = Trixi.BoundaryConditionDirichlet(bottom_boundary_value),)

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = initial_condition,
                             boundary_conditions = boundary_conditions,
                             domain = ((0.0,), (0.4,)), tspan = tspan)
end
