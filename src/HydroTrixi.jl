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
include("callbacks_step/amr/amr.jl")
include("analysis/analysis.jl")

include("equations/problems/problems.jl")
include("time_integration/time_integration.jl")

export HydrologicProblem
export RichardsEquation1D
export Haverkamp
export VanGenuchten
export water_content, pressure_head
export HydrologicProblemCelia1990
export HydrologicProblemRichardsClosedColumn
export BoundaryConditionDirichletPenalty
export SemidiscretizationImplicit
export TemporalOperatorConstitutive
export NoPassiveVariables
export PassiveVariables
export PassiveVariablesBoundaryFlux1D
export passive_variable_view
export passive_variables
export boundary_flux_integrals
export pressure_head_from_water_content
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
    plot_convergence_1d(ndofs, l2_errors, linf_errors; kwargs...)

Plot one-dimensional convergence data against the number of degrees of freedom, save the
figure to `output_path`, and return the `CairoMakie.Figure`.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_convergence_1d end

@doc raw"""
    animate_solution_1d(sol; kwargs...)
    animate_solution_1d(ode::SciMLBase.ODEProblem; callback, kwargs...)

Animate a one-dimensional solution and save it to `output_path`.

The first method renders frames from saved states in `sol`. The second method advances
an `ODEProblem` to the requested frame times while applying `callback`, which is useful
for animations with adaptive meshes.

These methods are provided by `HydroTrixiVisualizationExt` and become available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function animate_solution_1d end

function set_serif_tex_theme! end
function doubling_dof_ticks end
function plot_bottom_triangle! end
export plot_solution_1d
export animate_solution_1d
export plot_convergence_1d

end # module HydroTrixi
