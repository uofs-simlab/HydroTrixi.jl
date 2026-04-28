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
using SciMLBase
using Trixi
using LinearAlgebra

include("auxiliary/auxiliary.jl")
include("equations/equations.jl")
include("solvers/solvers.jl")
include("semidiscretization/semidiscretization.jl")

include("equations/problems/problems.jl")
include("time_integration/time_integration.jl")
include("visualization/visualization.jl")

export LinearDiffusionEquation1D
export HydrologicProblem
export RichardsEquation1D
export Haverkamp
export VanGenuchten
export HydrologicProblemCelia1990
export BoundaryConditionDirichletPenalty
export SemidiscretizationConstitutive
export default_algorithm
export plot_solution_1d
export animate_solution_1d
export compute_eoc, plot_convergence_1d
export examples_dir

end # module DiffuSEM
