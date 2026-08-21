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
default_stepsize_controller
solve_implicit
state_variable_norm
AbstractJacobianStrategy
DenseJacobian
SparseJacobian
compute_eoc
```

## Equations and problem setup

```@docs
HydrologicProblem
RichardsEquation1D
water_content
water_content_timederivative
mass_bias
mass_bias_history
water_capacity
hydraulic_conductivity
pressure_head
pressure_head_from_water_content
Haverkamp
VanGenuchten
HydrologicProblemCelia1990
HydrologicProblemRichardsManufacturedSolution
richards_manufactured_solution
source_terms_richards_manufactured_solution
HydrologicProblemRichardsClosedColumn
BoundaryConditionDirichletPenalty
MixedForm
PressureHeadForm
SemidiscretizationImplicit
semidiscretize
AbstractTemporalOperator
TemporalOperatorStandard
TemporalOperatorConstitutive
TemporalOperatorCapacity
HydroTrixi.AbstractPassiveVariables
NoPassiveVariables
PassiveVariablesBoundaryFlux1D
```

## Visualization

```@docs
plot_solution_1d
animate_solution_1d
plot_convergence_1d
plot_mass_bias
```
