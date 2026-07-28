# Run and plot a fixed-mesh refinement study for the Richards MMS benchmark

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf, @sprintf
using Trixi

elixir = joinpath(@__DIR__, "elixir_richards_manufactured_solution.jl")

final_time = 120.0
base_initial_refinement_level = 4
iterations = 4
polydeg = 3
reltol = 1.0e-9
abstol = 1.0e-11
dt = 1.0e-2
adaptive = true
saveat = Float64[]

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
result_prefix = "richards_manufactured_solution_refinement"
plot_format = "pdf"
run_refinement_study = true
plot_results = true

schemes = ((; name = "pressure_head", form = PressureHeadForm(),
            label_l2 = LaTeXString("Pressure-head \$L^2\$"),
            label_linf = LaTeXString("Pressure-head \$L^\\infty\$"), color = 2,
            marker = :rect, markersize = 11.0),
           (; name = "mixed", form = MixedForm(),
            label_l2 = LaTeXString("Mixed \$L^2\$"),
            label_linf = LaTeXString("Mixed \$L^\\infty\$"), color = 1, marker = :circle,
            markersize = 7.0))

function solve_refinement_level(elixir, scheme, level; final_time, polydeg, reltol, abstol,
                                dt, adaptive, saveat)
    Trixi.trixi_include(@__MODULE__, elixir; final_time = final_time,
                        initial_refinement_level = level, polydeg = polydeg,
                        form = scheme.form, reltol = reltol, abstol = abstol, dt = dt,
                        adaptive = adaptive, saveat = saveat, run_simulation = true)

    local_sol = Base.invokelatest(getproperty, @__MODULE__, :sol)
    local_analysis_callback = Base.invokelatest(getproperty, @__MODULE__,
                                                :analysis_callback)
    l2_error, linf_error = Base.invokelatest(local_analysis_callback, local_sol)

    return (; final_time = local_sol.t[end], l2_error = only(l2_error),
            linf_error = only(linf_error), retcode = local_sol.retcode)
end

function run_richards_refinement_study(elixir, schemes; base_initial_refinement_level,
                                       iterations, final_time, polydeg, reltol, abstol, dt,
                                       adaptive, saveat)
    levels = collect(base_initial_refinement_level:
                     (base_initial_refinement_level + iterations - 1))
    results = Dict{String, Any}()

    for scheme in schemes
        ndofs = Int[]
        final_times = Float64[]
        l2_errors = Float64[]
        linf_errors = Float64[]

        for level in levels
            println("Running Richards MMS $(scheme.name) refinement level $(level)")
            run_data = solve_refinement_level(elixir, scheme, level;
                                              final_time = final_time, polydeg = polydeg,
                                              reltol = reltol, abstol = abstol, dt = dt,
                                              adaptive = adaptive, saveat = saveat)
            push!(ndofs, (polydeg + 1) * 2^level)
            push!(final_times, run_data.final_time)
            push!(l2_errors, run_data.l2_error)
            push!(linf_errors, run_data.linf_error)

            println("$(scheme.name): retcode = $(run_data.retcode), " *
                    "ndofs = $(ndofs[end]), L2 = $(l2_errors[end]), " *
                    "Linf = $(linf_errors[end])")
        end

        results[scheme.name] = (; levels, ndofs, final_times, l2_errors, linf_errors)
    end

    return results
end

function write_refinement_table(output_path, schemes, results)
    open(output_path, "w") do io
        println(io,
                "form level ndofs final_time l2_pressure_head_error " *
                "linf_pressure_head_error l2_order linf_order")

        for scheme in schemes
            data = results[scheme.name]
            l2_orders = compute_eoc(data.l2_errors)
            linf_orders = compute_eoc(data.linf_errors)

            for i in eachindex(data.levels)
                l2_order = isfinite(l2_orders[i]) ? @sprintf("%.17e", l2_orders[i]) :
                           "missing"
                linf_order = isfinite(linf_orders[i]) ?
                             @sprintf("%.17e", linf_orders[i]) : "missing"
                @printf(io, "%s %d %d %.17e %.17e %.17e %s %s\n", scheme.name,
                        data.levels[i], data.ndofs[i], data.final_times[i],
                        data.l2_errors[i], data.linf_errors[i], l2_order, linf_order)
            end
        end
    end

    return output_path
end

function read_refinement_table(input_path, schemes)
    results = Dict{String, Any}()
    for scheme in schemes
        results[scheme.name] = (; levels = Int[], ndofs = Int[], final_times = Float64[],
                                l2_errors = Float64[], linf_errors = Float64[])
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
            form_name = values[1]
            haskey(results, form_name) || continue

            data = results[form_name]
            push!(data.levels, parse(Int, values[2]))
            push!(data.ndofs, parse(Int, values[3]))
            push!(data.final_times, parse(Float64, values[4]))
            push!(data.l2_errors, parse(Float64, values[5]))
            push!(data.linf_errors, parse(Float64, values[6]))
        end
    end

    return results
end

function plot_richards_refinement_study(schemes, results, output_path; polydeg,
                                        output_png_path = nothing)
    HydroTrixi.set_serif_tex_theme!()
    plot_font = HydroTrixi.DEFAULT_PLOT_FONT

    all_ndofs = Int[]
    for scheme in schemes
        append!(all_ndofs, results[scheme.name].ndofs)
    end
    fig = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
    xticks = HydroTrixi.doubling_dof_ticks(all_ndofs; base = minimum(all_ndofs))
    ax = Axis(fig[1, 1]; xlabel = LaTeXString("Degrees of freedom"),
              ylabel = LaTeXString("Pressure head error"), xlabelfont = plot_font,
              ylabelfont = plot_font, xticklabelfont = plot_font,
              yticklabelfont = plot_font, xscale = log10, yscale = log10, xticks = xticks)

    ax.xminorgridvisible = false
    ax.xminorticksvisible = false

    colors = Makie.wong_colors()
    for scheme in schemes
        data = results[scheme.name]
        scatterlines!(ax, data.ndofs, data.l2_errors; label = scheme.label_l2,
                      color = colors[scheme.color], linestyle = :solid, linewidth = 1.8,
                      marker = scheme.marker, markersize = scheme.markersize)
        scatterlines!(ax, data.ndofs, data.linf_errors; label = scheme.label_linf,
                      color = colors[scheme.color], linestyle = :dash, linewidth = 1.8,
                      marker = scheme.marker, markersize = scheme.markersize)
    end

    reference_data = results[first(schemes).name]
    y_ref = Inf
    for scheme in schemes
        data = results[scheme.name]
        y_ref = min(y_ref, data.l2_errors[end - 1], data.linf_errors[end - 1])
    end
    HydroTrixi.plot_bottom_triangle!(ax, reference_data.ndofs[end - 1],
                                     reference_data.ndofs[end], y_ref, polydeg + 1;
                                     triangle_shift = 1.5, trianglefontsize = 15,
                                     font = plot_font)

    axislegend(ax; position = (:left, :bottom), labelsize = 14, font = plot_font)

    mkpath(dirname(output_path))
    save(output_path, fig; px_per_unit = 1)
    if !isnothing(output_png_path) && output_png_path != output_path
        save(output_png_path, fig; px_per_unit = 1)
    end

    return fig
end

analysis_path = joinpath(plots_dir, "$(result_prefix)_analysis.dat")

if run_refinement_study
    results = run_richards_refinement_study(elixir, schemes;
                                            base_initial_refinement_level =
                                            base_initial_refinement_level,
                                            iterations = iterations,
                                            final_time = final_time, polydeg = polydeg,
                                            reltol = reltol, abstol = abstol, dt = dt,
                                            adaptive = adaptive, saveat = saveat)
    write_refinement_table(analysis_path, schemes, results)
    println("Saved Richards MMS refinement table to: $(analysis_path)")
else
    results = read_refinement_table(analysis_path, schemes)
end

if plot_results
    convergence_plot = joinpath(plots_dir, "$(result_prefix)_convergence.$(plot_format)")
    convergence_plot_png = joinpath(plots_dir, "$(result_prefix)_convergence.png")
    plot_richards_refinement_study(schemes, results, convergence_plot; polydeg = polydeg,
                                   output_png_path = convergence_plot_png)
    println("Saved Richards MMS refinement plot to: $(convergence_plot)")
    if convergence_plot_png != convergence_plot
        println("Saved Richards MMS refinement plot to: $(convergence_plot_png)")
    end
end
