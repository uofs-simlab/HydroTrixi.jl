"""
    HydroTrixi

**HydroTrixi.jl** is an adaptive discontinuous spectral-element framework for hydrologic
problems based on the Trixi.jl and SciML ecosystems.
"""
module HydroTrixi

using MuladdMacro
using OrdinaryDiffEq
using SciMLBase
using Trixi
using LinearAlgebra

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization.jl")
include("analysis/analysis.jl")

include("equations/problems/problems.jl")
include("time_integration/time_integration.jl")

export LinearDiffusionEquation1D
export HydrologicProblem
export RichardsEquation1D
export Haverkamp
export VanGenuchten
export HydrologicProblemCelia1990
export HydrologicProblemRichardsClosedColumn
export BoundaryConditionDirichletPenalty
export SemidiscretizationImplicit
export TemporalOperatorConstitutive
export default_algorithm
export compute_eoc
export examples_dir

# Visualization methods are added by HydroTrixiVisualizationExt when CairoMakie is
# loaded, to avoid making CairoMakie a hard dependency
const DEFAULT_PLOT_FONT = "CMU Serif"
const DEFAULT_SOLUTION_FIGSIZE = (500, 350)
const DEFAULT_CONVERGENCE_FIGSIZE = (500, 350)
function plot_solution_1d end
function animate_solution_1d end
function plot_convergence_1d end
function set_serif_tex_theme! end
function doubling_dof_ticks end
function plot_bottom_triangle! end
export plot_solution_1d
export animate_solution_1d
export plot_convergence_1d

end # module HydroTrixi
