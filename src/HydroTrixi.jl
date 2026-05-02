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
@doc raw"""
    plot_solution_1d(sol; kwargs...)

Plot the final one-dimensional solution profile stored in `sol`, save it to
`output_path`, and return the `CairoMakie.Figure`.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_solution_1d end

@doc raw"""
    animate_solution_1d(sol; kwargs...)

Animate the one-dimensional solution history stored in `sol`, save it to `output_path`,
and return the output path.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function animate_solution_1d end

@doc raw"""
    plot_convergence_1d(ndofs, l2_errors, linf_errors; kwargs...)

Plot one-dimensional convergence data against the number of degrees of freedom, save the
figure to `output_path`, and return the `CairoMakie.Figure`.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_convergence_1d end
function set_serif_tex_theme! end
function doubling_dof_ticks end
function plot_bottom_triangle! end
export plot_solution_1d
export animate_solution_1d
export plot_convergence_1d

end # module HydroTrixi
