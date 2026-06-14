# Run and plot a fixed-mesh time-step study for the Richards MMS benchmark

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf, @sprintf
using Trixi

elixir = joinpath(@__DIR__, "elixir_richards_manufactured_solution.jl")

final_time = 120.0
initial_refinement_level = 8
polydeg = 3
initial_dt = 4.0
iterations = 4
time_steps = initial_dt ./ 2.0 .^ collect(0:(iterations - 1))
reltol = 1.0e-9
abstol = 1.0e-11
saveat = Float64[]

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
result_prefix = "richards_manufactured_solution_temporal_convergence"
plot_format = "pdf"
run_temporal_study = true
plot_results = true
show_run_output = false
plot_ylims = (1.0e-12, 1.0e-3)
triangle_side = "right"
triangle_shift = 1.5

schemes = ((; name = "pressure_head",
            form = PressureHeadForm(),
            label_l2 = LaTeXString("Pressure-head \$L^2\$"),
            label_linf = LaTeXString("Pressure-head \$L^\\infty\$"),
            color = 2,
            marker = :rect,
            markersize = 11.0),
           (; name = "mixed",
            form = MixedForm(),
            label_l2 = LaTeXString("Mixed \$L^2\$"),
            label_linf = LaTeXString("Mixed \$L^\\infty\$"),
            color = 1,
            marker = :circle,
            markersize = 7.0))

function compute_time_eoc(errors, time_steps)
    length(errors) == length(time_steps) ||
        throw(ArgumentError("`errors` and `time_steps` must have equal length."))

    values = similar(errors, Float64)
    fill!(values, NaN)
    for i in (firstindex(errors) + 1):lastindex(errors)
        values[i] = log(errors[i - 1] / errors[i]) /
                    log(time_steps[i - 1] / time_steps[i])
    end

    return values
end

function run_time_step(elixir, scheme, time_step;
                       final_time,
                       initial_refinement_level,
                       polydeg,
                       reltol,
                       abstol,
                       saveat,
                       show_run_output)
    run_include = () -> Trixi.trixi_include(@__MODULE__, elixir;
                                            final_time = final_time,
                                            initial_refinement_level =
                                            initial_refinement_level,
                                            polydeg = polydeg,
                                            form = scheme.form,
                                            dt = time_step,
                                            adaptive = false,
                                            reltol = reltol,
                                            abstol = abstol,
                                            saveat = saveat,
                                            run_simulation = true)

    if show_run_output
        run_include()
    else
        redirect_stdout(devnull) do
            run_include()
        end
    end

    local_sol = Base.invokelatest(getproperty, @__MODULE__, :sol)
    local_analysis_callback = Base.invokelatest(getproperty, @__MODULE__,
                                                :analysis_callback)
    l2_error, linf_error = Base.invokelatest(local_analysis_callback, local_sol)

    return (; final_time = local_sol.t[end],
            nsteps = local_sol.destats.naccept,
            l2_error = only(l2_error),
            linf_error = only(linf_error),
            retcode = local_sol.retcode)
end

function run_richards_temporal_study(elixir, schemes, time_steps;
                                     final_time,
                                     initial_refinement_level,
                                     polydeg,
                                     reltol,
                                     abstol,
                                     saveat,
                                     show_run_output)
    results = Dict{String, Any}()

    for scheme in schemes
        final_times = Float64[]
        nsteps = Int[]
        l2_errors = Float64[]
        linf_errors = Float64[]
        retcodes = String[]

        for time_step in time_steps
            println("Running Richards MMS $(scheme.name) with dt = $(time_step)")
            run_data = run_time_step(elixir, scheme, time_step;
                                     final_time = final_time,
                                     initial_refinement_level =
                                     initial_refinement_level,
                                     polydeg = polydeg,
                                     reltol = reltol,
                                     abstol = abstol,
                                     saveat = saveat,
                                     show_run_output = show_run_output)
            push!(final_times, run_data.final_time)
            push!(nsteps, run_data.nsteps)
            push!(l2_errors, run_data.l2_error)
            push!(linf_errors, run_data.linf_error)
            push!(retcodes, string(run_data.retcode))

            println("$(scheme.name): retcode = $(run_data.retcode), " *
                    "steps = $(run_data.nsteps), L2 = $(l2_errors[end]), " *
                    "Linf = $(linf_errors[end])")
        end

        results[scheme.name] = (; time_steps = collect(time_steps),
                                final_times,
                                nsteps,
                                l2_errors,
                                linf_errors,
                                retcodes)
    end

    return results
end

function write_temporal_table(output_path, schemes, results)
    open(output_path, "w") do io
        println(io,
                "form dt nsteps final_time l2_pressure_head_error " *
                "linf_pressure_head_error l2_order linf_order retcode")

        for scheme in schemes
            data = results[scheme.name]
            l2_orders = compute_time_eoc(data.l2_errors, data.time_steps)
            linf_orders = compute_time_eoc(data.linf_errors, data.time_steps)

            for i in eachindex(data.time_steps)
                l2_order = isfinite(l2_orders[i]) ? @sprintf("%.17e", l2_orders[i]) :
                           "missing"
                linf_order = isfinite(linf_orders[i]) ?
                             @sprintf("%.17e", linf_orders[i]) : "missing"
                @printf(io,
                        "%s %.17e %d %.17e %.17e %.17e %s %s %s\n",
                        scheme.name,
                        data.time_steps[i],
                        data.nsteps[i],
                        data.final_times[i],
                        data.l2_errors[i],
                        data.linf_errors[i],
                        l2_order,
                        linf_order,
                        data.retcodes[i])
            end
        end
    end

    return output_path
end

function read_temporal_table(input_path, schemes)
    isfile(input_path) ||
        throw(ArgumentError("`input_path` must refer to an existing file."))

    results = Dict{String, Any}()
    for scheme in schemes
        results[scheme.name] = (; time_steps = Float64[],
                                final_times = Float64[],
                                nsteps = Int[],
                                l2_errors = Float64[],
                                linf_errors = Float64[],
                                retcodes = String[])
    end

    open(input_path, "r") do io
        header_read = false
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue

            if !header_read
                header_read = true
                continue
            end

            values = split(stripped_line)
            length(values) >= 9 ||
                throw(ArgumentError("`input_path` has a row with too few columns."))
            form_name = values[1]
            haskey(results, form_name) || continue

            data = results[form_name]
            push!(data.time_steps, parse(Float64, values[2]))
            push!(data.nsteps, parse(Int, values[3]))
            push!(data.final_times, parse(Float64, values[4]))
            push!(data.l2_errors, parse(Float64, values[5]))
            push!(data.linf_errors, parse(Float64, values[6]))
            push!(data.retcodes, values[9])
        end
    end

    return results
end

function plot_positive_order_triangle!(ax, x_left, x_right, y_ref, order;
                                       triangle_shift = 1.5,
                                       trianglefontsize = 12,
                                       font = HydroTrixi.DEFAULT_PLOT_FONT)
    x_right > x_left || throw(ArgumentError("`x_right` must exceed `x_left`."))
    ratio = x_right / x_left
    y_left = y_ref / triangle_shift
    y_right = y_left * ratio^order

    lines!(ax,
           [x_left, x_right, x_right, x_left],
           [y_left, y_left, y_right, y_left];
           color = :black)

    xm = 10^((log10(x_left) + 2 * log10(x_right)) / 3)
    ym = 10^((2 * log10(y_left) + log10(y_right)) / 3)
    text!(ax,
          xm,
          ym;
          text = string(order, ":1"),
          align = (:center, :center),
          color = :black,
          fontsize = trianglefontsize,
          font = font)

    return nothing
end

function time_step_ticks(values)
    ticks = Float64.(sort(unique(values)))
    labels = [@sprintf("%.4g", tick) for tick in ticks]
    return ticks, labels
end

function triangle_time_step_pair(sorted_time_steps, triangle_side)
    length(sorted_time_steps) >= 2 ||
        throw(ArgumentError("At least two time steps are required."))

    if triangle_side == "left"
        return sorted_time_steps[1], sorted_time_steps[2]
    elseif triangle_side == "right"
        return sorted_time_steps[end - 1], sorted_time_steps[end]
    end

    throw(ArgumentError("`triangle_side` must be \"left\" or \"right\"."))
end

function plot_richards_temporal_study(schemes, results, output_path;
                                      triangle_order = 5,
                                      triangle_side = "right",
                                      triangle_shift = 1.5,
                                      ylims = nothing,
                                      output_png_path = nothing)
    HydroTrixi.set_serif_tex_theme!()
    plot_font = HydroTrixi.DEFAULT_PLOT_FONT

    all_time_steps = Float64[]
    for scheme in schemes
        append!(all_time_steps, results[scheme.name].time_steps)
    end
    isempty(all_time_steps) && throw(ArgumentError("No temporal data found."))

    fig = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
    ax = Axis(fig[1, 1];
              xlabel = L"$\Delta t$ (s)",
              ylabel = LaTeXString("Pressure head error"),
              xlabelfont = plot_font,
              ylabelfont = plot_font,
              xticklabelfont = plot_font,
              yticklabelfont = plot_font,
              xscale = log10,
              yscale = log10,
              xticks = time_step_ticks(all_time_steps))

    ax.xminorgridvisible = false
    ax.xminorticksvisible = false
    isnothing(ylims) || ylims!(ax, ylims)

    colors = Makie.wong_colors()
    for scheme in schemes
        data = results[scheme.name]
        order = sortperm(data.time_steps)
        scatterlines!(ax,
                      data.time_steps[order],
                      data.l2_errors[order];
                      label = scheme.label_l2,
                      color = colors[scheme.color],
                      linestyle = :solid,
                      linewidth = 1.8,
                      marker = scheme.marker,
                      markersize = scheme.markersize)
        scatterlines!(ax,
                      data.time_steps[order],
                      data.linf_errors[order];
                      label = scheme.label_linf,
                      color = colors[scheme.color],
                      linestyle = :dash,
                      linewidth = 1.8,
                      marker = scheme.marker,
                      markersize = scheme.markersize)
    end

    minimum(length(results[scheme.name].time_steps) for scheme in schemes) >= 2 ||
        throw(ArgumentError("At least two time steps are required."))

    sorted_time_steps = sort(unique(all_time_steps))
    x_left, x_right = triangle_time_step_pair(sorted_time_steps,
                                              triangle_side)
    y_ref = Inf
    for scheme in schemes
        data = results[scheme.name]
        i_left = findfirst(==(x_left), data.time_steps)
        isnothing(i_left) && continue
        for error in (data.l2_errors[i_left], data.linf_errors[i_left])
            isfinite(error) && (y_ref = min(y_ref, error))
        end
    end
    isfinite(y_ref) ||
        throw(ArgumentError("No finite error value found for the convergence triangle."))
    plot_positive_order_triangle!(ax,
                                  x_left,
                                  x_right,
                                  y_ref,
                                  triangle_order;
                                  triangle_shift = triangle_shift,
                                  trianglefontsize = 15,
                                  font = plot_font)

    axislegend(ax; position = (:right, :bottom), labelsize = 14, font = plot_font)

    mkpath(dirname(output_path))
    save(output_path, fig; px_per_unit = 1)
    if !isnothing(output_png_path) && output_png_path != output_path
        save(output_png_path, fig; px_per_unit = 1)
    end

    return fig
end

analysis_path = joinpath(plots_dir, "$(result_prefix)_analysis.dat")

if run_temporal_study
    results = run_richards_temporal_study(elixir, schemes, time_steps;
                                          final_time = final_time,
                                          initial_refinement_level =
                                          initial_refinement_level,
                                          polydeg = polydeg,
                                          reltol = reltol,
                                          abstol = abstol,
                                          saveat = saveat,
                                          show_run_output = show_run_output)
    write_temporal_table(analysis_path, schemes, results)
    println("Saved Richards MMS temporal table to: $(analysis_path)")
else
    results = read_temporal_table(analysis_path, schemes)
end

if plot_results
    convergence_plot = joinpath(plots_dir,
                                "$(result_prefix).$(plot_format)")
    convergence_plot_png = joinpath(plots_dir, "$(result_prefix).png")
    plot_richards_temporal_study(schemes, results, convergence_plot;
                                 triangle_side = triangle_side,
                                 triangle_shift = triangle_shift,
                                 ylims = plot_ylims,
                                 output_png_path = convergence_plot_png)
    println("Saved Richards MMS temporal plot to: $(convergence_plot)")
    if convergence_plot_png != convergence_plot
        println("Saved Richards MMS temporal plot to: $(convergence_plot_png)")
    end
end
