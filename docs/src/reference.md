```@meta
CurrentModule = HydroTrixi
EditURL = "https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/docs/src/reference.md"
```

# Reference

## Package

```@docs
HydroTrixi
examples_dir
default_algorithm
compute_eoc
```

## Equations and problem setup

```@docs
HydrologicProblem
RichardsEquation1D
water_content
pressure_head_from_water_content
Haverkamp
VanGenuchten
HydrologicProblemCelia1990
HydrologicProblemRichardsClosedColumn
BoundaryConditionDirichletPenalty
SemidiscretizationImplicit
TemporalOperatorConstitutive
NoPassiveVariables
PassiveVariables
PassiveVariablesBoundaryFlux1D
passive_variable_view
passive_variables
boundary_flux_integrals
```

## Visualization

```@docs
plot_solution_1d
animate_solution_1d
plot_convergence_1d
```
