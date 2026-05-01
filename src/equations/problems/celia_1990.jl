@doc raw"""
    HydrologicProblemCelia1990()

Return the one-dimensional Richards-equation infiltration benchmark introduced by Celia,
Bouloutas, and Zarba (1990) as a `HydrologicProblem` for HydroTrixi semidiscretizations.

The benchmark uses the Haverkamp constitutive laws in the form reported by Ireson et al.
(2023), Eq. (25):
```math
\theta(\psi) = \theta_r + (\theta_s - \theta_r)
\frac{\alpha}{\alpha + |\psi|^\beta},
\qquad
K(\psi) = K_s\,\frac{A}{A + |\psi|^\gamma},
```
for ``\psi < 0``, with saturated values ``\theta = \theta_s`` and ``K = K_s`` for
``\psi \ge 0``. The setup uses depth ``z`` in metres, positive downward on
``z \in [0, 0.4]``, and time in seconds on ``t \in [0, 360]``. The pressure head is
initialized at ``-0.615`` m, with a Dirichlet boundary value of ``-0.207`` m at the top
boundary (`x_neg`) and ``-0.615`` m at the bottom boundary (`x_pos`).

The returned problem setup contains the fields `equations`, `constitutive_relation`,
`initial_condition`, `boundary_conditions`, `domain`, and `tspan`.

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
function HydrologicProblemCelia1990()
    soil_model = Haverkamp(saturated_hydraulic_conductivity = 9.44e-5,
                           alpha = 0.01936848004,
                           beta = 3.96,
                           A = 3.890790677e-4,
                           gamma = 4.74,
                           theta_s = 0.287,
                           theta_r = 0.075)
    equations = RichardsEquation1D(soil_model = soil_model)
    constitutive_relation(psi, equations) = water_content(psi, equations)
    initial_condition(x, t, equations) = Trixi.SVector(-0.615)
    top_boundary_value(x, t, equations) = Trixi.SVector(-0.207)
    bottom_boundary_value(x, t, equations) = Trixi.SVector(-0.615)
    boundary_conditions = (;
                           x_neg = Trixi.BoundaryConditionDirichlet(top_boundary_value),
                           x_pos = Trixi.BoundaryConditionDirichlet(bottom_boundary_value),)

    return HydrologicProblem(equations = equations,
                             constitutive_relation = constitutive_relation,
                             initial_condition = initial_condition,
                             boundary_conditions = boundary_conditions,
                             domain = ((0.0,), (0.4,)),
                             tspan = (0.0, 360.0))
end
