@muladd begin
#! format: noindent

@doc raw"""
    HydrologicProblemCeliaNewMexico(; tspan = (0.0, 86_400.0), penalty_factor = 1)

Return the second one-dimensional Richards-equation infiltration problem considered by
Celia, Bouloutas, and Zarba (1990), based on field data from New Mexico.

The soil follows the Mualem-van Genuchten constitutive laws
```math
S_{\mathrm{e}}(\psi) = \left(1 + (\alpha |\psi|)^n\right)^{-m},
\qquad
\kappa(\psi) = \kappa_{\mathrm{s}} S_{\mathrm{e}}^l
\left(1 - \left(1 - S_{\mathrm{e}}^{1/m}\right)^m\right)^2,
```
where ``m = 1 - 1/n`` and
``\vartheta = \theta_{\mathrm{r}} +
(\theta_{\mathrm{s}} - \theta_{\mathrm{r}})S_{\mathrm{e}}``. The parameter values are
``\alpha = 3.35\,\mathrm{m}^{-1}``, ``n = 2``, ``l = 0.5``,
``\kappa_{\mathrm{s}} = 9.22 \times 10^{-5}\,\mathrm{m}\,\mathrm{s}^{-1}``,
``\theta_{\mathrm{s}} = 0.368``, and ``\theta_{\mathrm{r}} = 0.102``.

Depth ``z`` is measured in metres, positive downward on ``z \in [0, 1]``, and time is
measured in seconds on ``t \in [0, 86400]`` (one day). The initial pressure head is
``-10`` m. Fixed pressure heads of ``-0.75`` m and ``-10`` m are imposed at the soil
surface (`x_neg`) and column bottom (`x_pos`), respectively. This one-metre domain follows
the depth axis used for the final profiles in Figure 3 of Celia et al. (1990); the paper's
nearby statement placing a boundary at 60 cm is inconsistent with that figure.

The Dirichlet boundaries use [`BoundaryConditionDirichletPenalty`](@ref). The dimensionless
`penalty_factor` is the coefficient ``C_\tau`` in the boundary penalty; setting it to zero
omits the additional divergence-flux penalty.

# References
- Celia, M. A., Bouloutas, E. T., Zarba, R. L. (1990). A general
  mass-conservative numerical solution for the unsaturated flow equation.
  *Water Resources Research*, 26(7), 1483-1496.
  [DOI: 10.1029/WR026i007p01483](https://doi.org/10.1029/WR026i007p01483)
"""
function HydrologicProblemCeliaNewMexico(; tspan = (0.0, 86_400.0), penalty_factor = 1)
    constitutive_model = VanGenuchten(saturated_hydraulic_conductivity = 9.22e-5,
                              alpha = 3.35, n = 2.0, pore_connectivity = 0.5,
                              theta_s = 0.368, theta_r = 0.102)
    equations = RichardsEquation1D(constitutive_model = constitutive_model)
    state_to_evolved = water_content
    evolved_to_state = pressure_head_from_water_content
    initial_condition(x, t, equations) = Trixi.SVector(-10.0)
    top_boundary_value(x, t, equations) = Trixi.SVector(-0.75)
    bottom_boundary_value(x, t, equations) = Trixi.SVector(-10.0)
    boundary_conditions = (;
                           x_neg = BoundaryConditionDirichletPenalty(top_boundary_value;
                                                                     penalty_factor),
                           x_pos = BoundaryConditionDirichletPenalty(bottom_boundary_value;
                                                                     penalty_factor))

    return HydrologicProblem(equations = equations, state_to_evolved = state_to_evolved,
                             evolved_to_state = evolved_to_state,
                             initial_condition = initial_condition,
                             boundary_conditions = boundary_conditions,
                             domain = ((0.0,), (1.0,)), tspan = tspan)
end
end # @muladd
