"""
    HydroTrixi

**HydroTrixi.jl** is an adaptive discontinuous spectral-element framework for hydrologic
problems based on the Trixi.jl and SciML ecosystems. Its one-dimensional Richards solver
supports mixed and pressure-head formulations with local discontinuous Galerkin spatial
discretization, adaptive implicit time integration, and adaptive mesh refinement.
"""
module HydroTrixi

import DiffEqBase
using MuladdMacro
import OrdinaryDiffEqCore
import OrdinaryDiffEqNonlinearSolve
import OrdinaryDiffEqRosenbrock
using PreallocationTools: GeneralLazyBufferCache
using SciMLBase
using Trixi
using LinearAlgebra
using LinearSolve
using SparseArrays

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization.jl")
include("callbacks_step/amr/amr.jl")
include("analysis/analysis.jl")
include("visualization/visualization.jl")

include("equations/problems/problems.jl")
include("time_integration/time_integration.jl")
include("time_integration/stepsize_controller_mapped_error.jl")

export HydrologicProblem
export RichardsEquation1D
export Haverkamp
export VanGenuchten
export AnalysisCallbackFullState
export water_content, water_content_timederivative
export mass_balance, mass_bias, mass_bias_history
export effective_saturation, water_capacity, hydraulic_conductivity, pressure_head
export HydrologicProblemCeliaHaverkamp
export HydrologicProblemCeliaNewMexico
export HydrologicProblemRichardsManufacturedSolution
export HydrologicProblemRichardsClosedColumn
export BoundaryConditionDirichletPenalty
export MixedForm
export PressureHeadForm
export AbstractJacobianStrategy
export DenseJacobian
export SparseJacobian
export SemidiscretizationImplicit
export AbstractTemporalOperator
export TemporalOperatorStandard
export TemporalOperatorConstitutive
export TemporalOperatorCapacity
export evolved_variable_block
export state_variable_block
export AbstractPassiveVariables
export NoPassiveVariables
export PassiveVariablesBoundaryFlux1D
export pressure_head_from_water_content
export default_algorithm
export default_stepsize_controller
export solve_implicit
export pressure_head_out_of_domain
export ControllerTwoThreshold
export IndicatorTotalVariation
export compute_eoc
export examples_dir

end # module HydroTrixi
