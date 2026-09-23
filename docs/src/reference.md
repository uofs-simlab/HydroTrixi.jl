```@meta
CurrentModule = HydroTrixi
EditURL = "https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/docs/src/reference.md"
```

# Reference

## Package

```@docs
HydroTrixi
examples_dir
```

## Problem setup and boundary conditions

```@docs
HydrologicProblem
HydrologicProblemCeliaHaverkamp
HydrologicProblemCeliaNewMexico
HydrologicProblemRichardsManufacturedSolution
HydrologicProblemRichardsClosedColumn
BoundaryConditionDirichletPenalty
```

## Equations and constitutive models

```@docs
RichardsEquation1D
Haverkamp
VanGenuchten
water_content
water_content_timederivative
effective_saturation
water_capacity
hydraulic_conductivity
pressure_head
pressure_head_from_water_content
```

## Semi-discrete formulations and temporal operators

```@docs
SemidiscretizationImplicit
semidiscretize
MixedForm
PressureHeadForm
AbstractTemporalOperator
TemporalOperatorStandard
TemporalOperatorConstitutive
TemporalOperatorCapacity
evolved_variable_block
state_variable_block
```

### Passive diagnostic variables

```@docs
AbstractPassiveVariables
NoPassiveVariables
PassiveVariablesBoundaryFlux1D
```

### Jacobian strategies

```@docs
AbstractJacobianStrategy
DenseJacobian
SparseJacobian
```

## Spatial adaptivity

```@docs
ControllerTwoThreshold
IndicatorTotalVariation
```

## Time integration

```@docs
default_algorithm
default_stepsize_controller
solve_implicit
pressure_head_out_of_domain
```

## Analysis

```@docs
AnalysisCallbackFullState
mass_balance
mass_bias
mass_bias_history
compute_eoc
```

## Visualization

```@docs
plot_solution_1d
animate_solution_1d
plot_convergence_1d
plot_mass_bias
```
