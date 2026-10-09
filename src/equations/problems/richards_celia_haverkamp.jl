@muladd begin
#! format: noindent

@doc raw"""
    HydrologicProblemCeliaHaverkamp(; tspan = (0.0, 360.0), penalty_factor = 1)

Return the Haverkamp infiltration problem,
the first one-dimensional Richards equation infiltration problem
introduced by Celia, Bouloutas, and Zarba (1990),
as a `HydrologicProblem` for HydroTrixi.jl semidiscretizations.

The problem uses the constitutive model of Haverkamp et al. (1977), given by
```math
\vartheta(\psi) = \theta_{\mathrm{r}} +
                  (\theta_{\mathrm{s}} - \theta_{\mathrm{r}})
\frac{1}{1 + (a_{\mathrm{H}}|\psi|)^{\beta_{\mathrm{H}}}},
\qquad
\kappa(\psi) = \frac{\kappa_{\mathrm{s}}}
{1 + (b_{\mathrm{H}}|\psi|)^{\gamma_{\mathrm{H}}}},
```
for ``\psi < 0``, with saturated values ``\vartheta(\psi) = \theta_{\mathrm{s}}``
and ``\kappa(\psi) = \kappa_{\mathrm{s}}`` for ``\psi \ge 0``.
The parameters ``a_{\mathrm{H}}``, ``b_{\mathrm{H}}``, ``\beta_{\mathrm{H}}``,
and ``\gamma_{\mathrm{H}}`` correspond to the [`Haverkamp`](@ref) keywords
`a`, `b`, `beta`, and `gamma`, respectively.

The setup uses depth ``z`` in metres, positive downward on ``z \in [0, 0.4]``,
and time in seconds, with the default interval ``t \in [0, 360]``.
The pressure head is initialized as ``\psi(z,0)=-0.615`` m.
The Dirichlet data are ``\psi_{\mathrm{T}}(t)=-0.207`` m at the soil surface (`x_neg`)
and ``\psi_{\mathrm{B}}(t)=-0.615`` m at the bottom of the column (`x_pos`).

The Dirichlet boundaries use [`BoundaryConditionDirichletPenalty`](@ref).
The dimensionless `penalty_factor` is the coefficient ``C_\tau`` in the boundary penalty.
Its default value is one;
setting it to zero omits the additional divergence-flux penalty.

The returned problem setup contains the fields `equations`, `state_to_evolved`,
`evolved_to_state`, `initial_condition`, `boundary_conditions`, `domain`, and `tspan`.

# References
- Haverkamp, R., Vauclin, M., Touma, J., Wierenga, P. J., Vachaud, G. (1977).
  A comparison of numerical simulation models for one-dimensional infiltration.
  *Soil Science Society of America Journal*, 41(2), 285-294.
  [DOI: 10.2136/sssaj1977.03615995004100020024x](https://doi.org/10.2136/sssaj1977.03615995004100020024x)
- Celia, M. A., Bouloutas, E. T., Zarba, R. L. (1990).
  A general mass-conservative numerical solution for the unsaturated flow equation.
  *Water Resources Research*, 26(7), 1483-1496.
  [DOI: 10.1029/WR026i007p01483](https://doi.org/10.1029/WR026i007p01483)
"""
function HydrologicProblemCeliaHaverkamp(; tspan = (0.0, 360.0), penalty_factor = 1)
    constitutive_model = Haverkamp(saturated_hydraulic_conductivity = 9.44e-5,
                           a = 2.7073950541818448, beta = 3.96,
                           b = 5.2408447406427436, gamma = 4.74,
                           theta_s = 0.287, theta_r = 0.075)
    equations = RichardsEquation1D(constitutive_model = constitutive_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content
    initial_condition(x, t, equations) = Trixi.SVector(-0.615)
    top_boundary_value(x, t, equations) = Trixi.SVector(-0.207)
    bottom_boundary_value(x, t, equations) = Trixi.SVector(-0.615)
    boundary_conditions = (;
                           x_neg = BoundaryConditionDirichletPenalty(top_boundary_value;
                                                                     penalty_factor),
                           x_pos = BoundaryConditionDirichletPenalty(bottom_boundary_value;
                                                                     penalty_factor))

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = initial_condition,
                             boundary_conditions = boundary_conditions,
                             domain = ((0.0,), (0.4,)), tspan = tspan)
end
end # @muladd
