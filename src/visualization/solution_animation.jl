# Extract the curve and mesh vertices from the same PlotData1D snapshot
function animation_curve_data_1d(u_ode, semi; component)
    pd, local_component = plot_data_1d(u_ode, semi; component = component)
    mesh_vertices_x = isnothing(pd.mesh_vertices_x) ? Float64[] :
                      collect(pd.mesh_vertices_x)
    return collect(pd.x), plot_component_data(pd, local_component), mesh_vertices_x
end

# Store each curve as points so adaptive meshes update atomically
function solution_points_1d(x, y)
    return Makie.Point2f.(x, y)
end

# Create the shared Makie objects for one-dimensional animations
function initialize_animation_plot_1d(x, y, mesh_vertices_x, t; exact_solution,
                                      numerical_label, exact_label, xlabel, ylabel, font,
                                      size, fontsize, legendfontsize, linewidth, markersize,
                                      show_nodes, show_boundaries, xlabelfont, ylabelfont,
                                      titlefont, xticklabelfont, yticklabelfont, legendfont,
                                      legend_position, xlims, ylims)
    HydroTrixi.set_serif_tex_theme!(font = font)

    show_exact = !isnothing(exact_solution)
    points_obs = Observable(solution_points_1d(x, y))
    mesh_vertices_x_obs = Observable(mesh_vertices_x)

    fig = Figure(size = size, fontsize = fontsize)
    ax = solution_axis(fig; xlabel = xlabel, ylabel = ylabel, xlabelfont = xlabelfont,
                       ylabelfont = ylabelfont, titlefont = titlefont,
                       xticklabelfont = xticklabelfont, yticklabelfont = yticklabelfont,
                       xlims = xlims, ylims = ylims)

    if show_boundaries
        vlines!(ax, mesh_vertices_x_obs; color = (:gray, 0.45), linewidth = 0.9)
    end

    if show_nodes
        scatterlines!(ax, points_obs; label = series_label(show_exact, numerical_label),
                      linewidth = linewidth, markersize = markersize,
                      color = Makie.wong_colors()[1])
    else
        lines!(ax, points_obs; label = series_label(show_exact, numerical_label),
               linewidth = linewidth, color = Makie.wong_colors()[1])
    end

    x_exact = nothing
    y_exact_obs = nothing
    if show_exact
        x_exact = exact_solution_x(x)
        y_exact_obs = Observable(exact_solution_values(exact_solution, x_exact, t))
        lines!(ax, x_exact, y_exact_obs; label = exact_label, linewidth = linewidth,
               linestyle = :dash, color = Makie.wong_colors()[2],)
    end
    add_legend!(ax; position = legend_position, font = legendfont,
                labelsize = legendfontsize, show_legend = show_exact)

    return (; fig, points_obs, mesh_vertices_x_obs, x_exact, y_exact_obs)
end

# Update the animated observables for the next rendered frame
function update_animation_plot_1d!(animation_plot, x, y, mesh_vertices_x, t; exact_solution)
    animation_plot.points_obs[] = solution_points_1d(x, y)
    animation_plot.mesh_vertices_x_obs[] = mesh_vertices_x

    if !isnothing(animation_plot.y_exact_obs)
        animation_plot.y_exact_obs[] = exact_solution_values(exact_solution,
                                                             animation_plot.x_exact, t)
    end

    return nothing
end

# Create the animation destination directory when necessary
function prepare_animation_output_path(output_path)
    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)

    return nothing
end

@doc raw"""
    animate_solution_1d(sol; output_path = joinpath(pwd(), "solution_1d.mp4"),
                        component = 1, frame_indices = nothing, kwargs...)

Animate a one-dimensional solution history stored in `sol`, write the animation to
`output_path`, and return `output_path`.

The frames are taken from saved solution states. If `frame_indices` is `nothing`, every
saved state in `sol.t` is rendered; otherwise, only the selected indices are rendered.
Set `component` to choose the plotted variable. Optional keyword arguments control axis
labels and limits, figure size, fonts, line and marker styling, node markers, element
boundary guides, exact-solution overlays, and the output `framerate`.

`exact_solution`, when supplied, is called as `exact_solution(Trixi.SVector(x), t)`.
"""
function HydroTrixi.animate_solution_1d(sol;
                                        output_path = joinpath(pwd(), "solution_1d.mp4"),
                                        component = 1, exact_solution = nothing,
                                        numerical_label = LaTeXString("Numerical"),
                                        exact_label = LaTeXString("Exact"), xlabel = L"$x$",
                                        ylabel = L"$u(x,t)$",
                                        font = HydroTrixi.DEFAULT_PLOT_FONT,
                                        size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE,
                                        fontsize = 15, legendfontsize = 14, linewidth = 2.0,
                                        markersize = 7.0, show_nodes = false,
                                        show_element_boundaries = false,
                                        xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        titlefont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legendfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legend_position = :rb, xlims = nothing,
                                        ylims = nothing, framerate = 24,
                                        frame_indices = nothing,)
    indices = if isnothing(frame_indices)
        collect(eachindex(sol.t))
    else
        collect(frame_indices)
    end
    semi = sol.prob.p
    first_idx = first(indices)
    x, y, mesh_vertices_x = animation_curve_data_1d(sol.u[first_idx], semi;
                                                    component = component)
    animation_plot = initialize_animation_plot_1d(x, y, mesh_vertices_x, sol.t[first_idx];
                                                  exact_solution = exact_solution,
                                                  numerical_label = numerical_label,
                                                  exact_label = exact_label,
                                                  xlabel = xlabel, ylabel = ylabel,
                                                  font = font, size = size,
                                                  fontsize = fontsize,
                                                  legendfontsize = legendfontsize,
                                                  linewidth = linewidth,
                                                  markersize = markersize,
                                                  show_nodes = show_nodes,
                                                  show_boundaries = show_element_boundaries,
                                                  xlabelfont = xlabelfont,
                                                  ylabelfont = ylabelfont,
                                                  titlefont = titlefont,
                                                  xticklabelfont = xticklabelfont,
                                                  yticklabelfont = yticklabelfont,
                                                  legendfont = legendfont,
                                                  legend_position = legend_position,
                                                  xlims = xlims, ylims = ylims)

    prepare_animation_output_path(output_path)

    record(animation_plot.fig, output_path, indices; framerate = framerate,
           px_per_unit = 1) do i
        frame_data = animation_curve_data_1d(sol.u[i], semi; component = component)
        x_frame, y_frame, mesh_vertices_x_frame = frame_data
        update_animation_plot_1d!(animation_plot, x_frame, y_frame, mesh_vertices_x_frame,
                                  sol.t[i]; exact_solution = exact_solution)
    end

    return output_path
end

@doc raw"""
    animate_solution_1d(ode::SciMLBase.ODEProblem; callback,
                        output_path = joinpath(pwd(), "solution_1d.mp4"),
                        component = 1, saveat = range(first(ode.tspan),
                        last(ode.tspan); length = 121), kwargs...)

Animate a one-dimensional `ODEProblem` by advancing an integrator to each requested
frame time, write the animation to `output_path`, and return `output_path`.

This method is useful for adaptive meshes because the supplied `callback` is active while
frames are generated, so mesh adaptation and other step callbacks remain synchronized with
the plotted state. `saveat` gives the animation frame times and must be non-empty,
sorted, and contained in `ode.tspan`. `alg`, `dt`, `adaptive`, and additional SciML
integrator options control the integrator used for frame generation.
"""
function HydroTrixi.animate_solution_1d(ode::SciMLBase.ODEProblem; callback,
                                        output_path = joinpath(pwd(), "solution_1d.mp4"),
                                        component = 1,
                                        alg = HydroTrixi.default_algorithm(ode.p),
                                        dt = nothing, adaptive = true,
                                        saveat = range(first(ode.tspan), last(ode.tspan);
                                                       length = 121),
                                        exact_solution = nothing,
                                        numerical_label = LaTeXString("Numerical"),
                                        exact_label = LaTeXString("Exact"), xlabel = L"$x$",
                                        ylabel = L"$u(x,t)$",
                                        font = HydroTrixi.DEFAULT_PLOT_FONT,
                                        size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE,
                                        fontsize = 15, legendfontsize = 14, linewidth = 2.0,
                                        markersize = 7.0, show_nodes = false,
                                        show_element_boundaries = false,
                                        xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        titlefont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legendfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legend_position = :rb, xlims = nothing,
                                        ylims = nothing, framerate = 24, kwargs...)
    frame_times = collect(saveat)
    if !issorted(frame_times)
        throw(ArgumentError("`saveat` must be sorted."))
    end

    t_start, t_end = ode.tspan
    if first(frame_times) < t_start
        throw(ArgumentError("`saveat` starts before the ODE time span."))
    end
    if last(frame_times) > t_end
        throw(ArgumentError("`saveat` ends after the ODE time span."))
    end

    integrator = if isnothing(dt)
        SciMLBase.init(ode, alg; adaptive = adaptive, callback = callback,
                       Trixi.ode_default_options()..., kwargs...)
    else
        SciMLBase.init(ode, alg; dt = dt, adaptive = adaptive, callback = callback,
                       Trixi.ode_default_options()..., kwargs...)
    end

    for frame_time in frame_times
        if frame_time > integrator.t
            SciMLBase.add_tstop!(integrator, frame_time)
        end
    end

    semi = ode.p
    x, y, mesh_vertices_x = animation_curve_data_1d(integrator.u, semi;
                                                    component = component)
    animation_plot = initialize_animation_plot_1d(x, y, mesh_vertices_x, integrator.t;
                                                  exact_solution = exact_solution,
                                                  numerical_label = numerical_label,
                                                  exact_label = exact_label,
                                                  xlabel = xlabel, ylabel = ylabel,
                                                  font = font, size = size,
                                                  fontsize = fontsize,
                                                  legendfontsize = legendfontsize,
                                                  linewidth = linewidth,
                                                  markersize = markersize,
                                                  show_nodes = show_nodes,
                                                  show_boundaries = show_element_boundaries,
                                                  xlabelfont = xlabelfont,
                                                  ylabelfont = ylabelfont,
                                                  titlefont = titlefont,
                                                  xticklabelfont = xticklabelfont,
                                                  yticklabelfont = yticklabelfont,
                                                  legendfont = legendfont,
                                                  legend_position = legend_position,
                                                  xlims = xlims, ylims = ylims)

    prepare_animation_output_path(output_path)

    record(animation_plot.fig, output_path, frame_times; framerate = framerate,
           px_per_unit = 1) do frame_time
        tolerance = 100 * eps(max(abs(integrator.t), abs(frame_time), 1.0))
        while integrator.t < frame_time - tolerance
            SciMLBase.step!(integrator)
        end

        frame_data = animation_curve_data_1d(integrator.u, semi; component = component)
        x_frame, y_frame, mesh_vertices_x_frame = frame_data
        update_animation_plot_1d!(animation_plot, x_frame, y_frame, mesh_vertices_x_frame,
                                  integrator.t; exact_solution = exact_solution)
    end

    return output_path
end
