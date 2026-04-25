"""
    DiffuSEM

DiffuSEM.jl is a numerical solver package for parabolic partial differential equations 
based on the Trixi.jl and SciML ecosystems.
"""
module DiffuSEM

using CairoMakie
using LaTeXStrings
using MuladdMacro
using OrdinaryDiffEq
using Trixi

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("solvers/solvers.jl")
include("time_integration/time_integration.jl")
include("plotting/style.jl")
include("plotting/solution_plot.jl")
include("plotting/solution_animation.jl")
include("plotting/convergence_plot.jl")

export LinearDiffusionEquation1D
export AbstractDiffusionEquation1D
export BoundaryConditionDirichletPenalty
export SemidiscretizationParabolic
export default_algorithm
export plot_solution_1d
export animate_solution_1d
export compute_eoc, plot_convergence_1d
export examples_dir

end # module DiffuSEM
