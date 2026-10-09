# Visualization methods are added by HydroTrixiVisualizationExt when CairoMakie is
# loaded, to avoid making CairoMakie a hard dependency.
const DEFAULT_PLOT_FONT = "CMU Serif"
const DEFAULT_SOLUTION_FIGSIZE = (500, 350)
const DEFAULT_CONVERGENCE_FIGSIZE = (500, 350)

@doc raw"""
    plot_solution_1d(sol; kwargs...)
    plot_solution_1d(result::ImplicitSolveResult; index=lastindex(result.sol.u), kwargs...)
    plot_solution_1d(data::NamedTuple; kwargs...)

Plot a one-dimensional solution profile, save it to `output_path`,
and return the `CairoMakie.Figure`.
With an `ImplicitSolveResult`, `index` selects a saved state;
plotting an earlier state requires recorded mesh history.
The horizontal limits default to the domain endpoints.
Set `show_element_boundaries = true` to draw the mesh guides used in animations.

A named-tuple input must contain `time`, `x`, `values`, and `mesh_vertices_x`, as returned
by [`solution_data_1d`](@ref). This also supports profiles read from saved numerical tables.
Use the same axis, styling, and exact-solution keywords as for solution inputs.

```julia
using HydroTrixi, CairoMakie, LaTeXStrings

data = solution_data_1d(result; component = 2)
plot_solution_1d(data; output_path = "pressure_head.pdf", show_element_boundaries = true)
plot_solution_1d((; time = 0.0, x = [0.0, 0.1, 0.2], values = [-0.6, -0.5, -0.4],
                    mesh_vertices_x = [0.0, 0.2]); output_path = "profile.pdf")
```

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_solution_1d end

@doc raw"""
    plot_convergence_1d(series; kwargs...)

Plot one-dimensional convergence data, save the figure to `output_path`,
and return the `CairoMakie.Figure`.

Each group in `series` must provide a shared `x` vector, an `errors` collection containing
one or more error vectors, and matching `labels`. It may also provide shared `color`,
`marker`, and `markersize` values. Integer colours select entries from Makie's Wong palette.
Pass `triangle_order` to infer a reference triangle for groups with shared `x` values.
By default, `triangle_gap_factor = 1.5` places the triangle below the curves by that factor.
Set `triangle_position = :above` to place it above them. Use `triangle_slope = :positive`
for mesh spacings or time steps that decrease with refinement;
the default `:negative` is for increasing degrees of freedom.
Placement uses both errors at the finest interval.
The x ticks are inferred for doubling degrees of freedom; pass `xticks = nothing` or
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

Plot the absolute bias error ``|\epsilon_{\mathrm{b}}(t^n)|`` defined by
[`mass_bias`](@ref), using a logarithmic axis by default.
Read the bias history from `sol`, a Trixi.jl analysis file, or multiple sources.
Pass `yscale = identity` for a linear axis that includes zero.
Zero magnitudes are omitted from the logarithmic plot.
Save the plot to `output_path` and return the `CairoMakie.Figure`.
For multiple sources, `labels` controls the legend
while `colors` and `linestyles` select each curve's style.
Pass `xticks`, `yticks`, `xlims`, and `ylims` to set axis ticks and limits.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_mass_bias_magnitude end

@doc raw"""
    mass_bias_magnitude_axis(sources; time_column = "time",
                             mass_balance_column = "mass_balance")

Return logarithmic `limits` and `ticks` for the absolute bias error across all `sources`.
Sources may be solutions or analysis files, as for [`plot_mass_bias_magnitude`](@ref).
Use the returned values as `ylims` and `yticks` to share one axis across figures.

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function mass_bias_magnitude_axis end

@doc raw"""
    plot_time_steps(source; output_path = "time_steps.pdf", kwargs...)

Plot accepted-step sizes against their end times, save the figure to `output_path`,
and return the `CairoMakie.Figure`.
Both axes are linear by default; set `xscale` or `yscale` to change them.
Pass `xticks`, `yticks`, `xlims`, and `ylims` to set axis ticks and limits,
and `xlabel` and `ylabel` to change axis labels.

Supply an analysis-file path, a named tuple containing `times` and `dts`,
or a tuple or vector of these sources.
File inputs use [`accepted_step_history`](@ref), omitting step zero;
pass `step_column`, `time_column`, and `dt_column` for alternate column names.
Array histories are plotted as supplied. A saved ODE solution is not a history input:
its saved times need not include every accepted step.

For multiple sources, `labels`, `colors`, and `linestyles` each take one entry per source.
The legend is omitted when `labels` is not supplied.
Set its position with `legend_position`.

```julia
using HydroTrixi, CairoMakie, LaTeXStrings

plot_time_steps("analysis.dat"; output_path = "steps.pdf")
history = accepted_step_history("analysis.dat")
plot_time_steps((; times = history.times, dts = history.dts);
                output_path = "steps_from_arrays.pdf")
plot_time_steps(["conservation_1.dat", "conservation_2.dat"];
                step_column = "accepted_step", time_column = "time_s", dt_column = "dt_s",
                labels = ["Run 1", "Run 2"], output_path = "comparison.pdf")
```

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_time_steps end

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
@doc raw"""
    plot_reference_triangle!(ax, series_groups, order; position = :below,
                             triangle_slope = :negative, gap_factor = 1.5,
                             trianglefontsize = 12, font = DEFAULT_PLOT_FONT)
    plot_reference_triangle!(ax, coarse_x, fine_x, reference_error, order; kwargs...)

Add a reference triangle to a logarithmic convergence axis and return `nothing`.
The series groups use the same `x` and `errors` fields as [`plot_convergence_1d`](@ref),
with shared `x` values.
Span the finest refinement interval and use both endpoint errors
to place the sloping edge below all curves (`position = :below`)
or above all curves (`position = :above`).
`gap_factor` sets the multiplicative separation from the curves.

Use `triangle_slope = :positive` for decreasing mesh spacings or time steps and
`:negative` for increasing degrees of freedom. The explicit-endpoint method takes
`reference_error` at the coarse endpoint before applying the gap factor.

```julia
fig = plot_convergence_1d(groups; output_path = "convergence.pdf", xticks = nothing)
plot_reference_triangle!(fig.content[1], groups, 4; triangle_slope = :positive)
plot_reference_triangle!(fig.content[1], groups, 5;
                         triangle_slope = :positive, position = :above)
save("convergence.pdf", fig)
```

This method is provided by `HydroTrixiVisualizationExt` and becomes available when
`CairoMakie` and `LaTeXStrings` are loaded.
"""
function plot_reference_triangle! end

export plot_solution_1d
export animate_solution_1d
export plot_convergence_1d
export plot_mass_bias_magnitude
export mass_bias_magnitude_axis
export plot_time_steps
export plot_reference_triangle!
