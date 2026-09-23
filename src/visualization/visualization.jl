# Visualization methods are added by HydroTrixiVisualizationExt when CairoMakie is
# loaded, to avoid making CairoMakie a hard dependency.
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
    plot_convergence_1d(series; kwargs...)

Plot one-dimensional convergence data, save the figure to `output_path`, and return the
`CairoMakie.Figure`.

Each group in `series` must provide a shared `x` vector, an `errors` collection containing
one or more error vectors, and matching `labels`. It may also provide shared `color`,
`marker`, and `markersize` values. Integer colors select entries from Makie's Wong palette.
Pass `triangle_order` to infer a reference triangle for groups with shared `x` values. By
default, `triangle_gap_factor = 1.5` places the triangle below the nearest curve by that
factor. The x ticks are inferred for doubling degrees of freedom; pass `xticks = nothing` or
explicit ticks for other x axes.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_convergence_1d end

@doc raw"""
    plot_mass_bias(sol; kwargs...)
    plot_mass_bias(analysis_path::AbstractString;
                   time_column = "time", mass_balance_column = "mass_balance", kwargs...)
    plot_mass_bias(sources; labels, kwargs...)

Plot the water mass bias obtained by differencing the saved history of
```math
F_K(t)-F_0(t)
-\sum_{k=1}^{K}J_k\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
 \boldsymbol{\theta}_k(t)
```
in `sol`, in a Trixi.jl analysis file, or in multiple sources. Save the plot to
`output_path` and return the `CairoMakie.Figure`. By default, the bias is scaled by an
automatically selected power of ten and the corresponding exponent is shown above the
top-left corner of the axis.
Pass `absolute = true, yscale = log10` to plot ``|\epsilon_{\mathrm{B}}|`` on a logarithmic
axis.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_mass_bias end

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
export plot_mass_bias
