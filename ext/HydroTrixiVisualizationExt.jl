module HydroTrixiVisualizationExt

using CairoMakie
using HydroTrixi
using LaTeXStrings
using SciMLBase
using Trixi

const DEFAULT_PLOT_FONT = HydroTrixi.DEFAULT_PLOT_FONT
const DEFAULT_SOLUTION_FIGSIZE = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE
const DEFAULT_CONVERGENCE_FIGSIZE = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE
const SemidiscretizationImplicit = HydroTrixi.SemidiscretizationImplicit
const evolved_variable_view = HydroTrixi.evolved_variable_view
const state_variable_view = HydroTrixi.state_variable_view
const compute_eoc = HydroTrixi.compute_eoc
const set_serif_tex_theme! = HydroTrixi.set_serif_tex_theme!
const doubling_dof_ticks = HydroTrixi.doubling_dof_ticks
const plot_bottom_triangle! = HydroTrixi.plot_bottom_triangle!

include("../src/visualization/visualization.jl")

end
