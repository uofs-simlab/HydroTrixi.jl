module HydroTrixiVisualizationExt

using CairoMakie
using HydroTrixi
using LaTeXStrings
using SciMLBase
using Trixi

include("../src/visualization/style.jl")
include("../src/visualization/solution_plot.jl")
include("../src/visualization/solution_animation.jl")
include("../src/visualization/convergence_plot.jl")
include("../src/visualization/mass_bias_plot.jl")

end
