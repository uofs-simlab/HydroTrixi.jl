# Visualization methods are added by HydroTrixiVisualizationExt when CairoMakie is
# loaded, to avoid making CairoMakie a hard dependency.
const DEFAULT_PLOT_FONT = "CMU Serif"
const DEFAULT_SOLUTION_FIGSIZE = (500, 350)
const DEFAULT_CONVERGENCE_FIGSIZE = (500, 350)

@doc raw"""
    plot_solution_1d(sol; kwargs...)
    plot_solution_1d(result::ImplicitSolveResult; index=lastindex(result.sol.u), kwargs...)

Plot a one-dimensional solution profile, save it to `output_path`, and return the
`CairoMakie.Figure`. With an `ImplicitSolveResult`, `index` selects a saved state;
plotting an earlier state requires recorded mesh history.
Set `show_element_boundaries = true` to draw the mesh guides used in animations.

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
explicit ticks for other x axes. Pass `yticks` to set the y-axis ticks.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_convergence_1d end

@doc raw"""
    plot_mass_bias_magnitude(sol; kwargs...)
    plot_mass_bias_magnitude(analysis_path::AbstractString;
                             time_column = "time", mass_balance_column = "mass_balance",
                             kwargs...)
    plot_mass_bias_magnitude(sources; labels, colors, linestyles, kwargs...)

Plot the absolute water mass bias on a logarithmic axis by default, obtained by
differencing the saved history of
```math
F_{\mathrm{B}}(t)-F_{\mathrm{T}}(t)
-\sum_{k=1}^{K}J_k\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
 \boldsymbol{\theta}_k(t)
```
in `sol`, in a Trixi.jl analysis file, or in multiple sources. Pass `yscale = identity`
for a linear axis that includes zero. Zero magnitudes are omitted from the logarithmic
plot. Save the plot to `output_path` and return the
`CairoMakie.Figure`. For multiple sources, `labels` controls the legend while `colors`
and `linestyles` select each curve's style. Pass `xticks`, `yticks`, `xlims`, and `ylims`
to set axis ticks and limits.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_mass_bias_magnitude end

@doc raw"""
    mass_bias_magnitude_axis(sources; time_column = "time",
                             mass_balance_column = "mass_balance")

Return logarithmic `limits` and `ticks` for the absolute mass bias across all `sources`.
Sources may be solutions or analysis files, as for [`plot_mass_bias_magnitude`](@ref).
Use the returned values as `ylims` and `yticks` to share one axis across figures.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function mass_bias_magnitude_axis end

@doc raw"""
    animate_solution_1d(sol; kwargs...)
    animate_solution_1d(result::ImplicitSolveResult; kwargs...)

Animate a saved one-dimensional solution and save it to `output_path`.
An `ImplicitSolveResult` supplies its recorded mesh history automatically.

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
export plot_mass_bias_magnitude
export mass_bias_magnitude_axis
