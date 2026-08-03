# Run and plot a time-step study for the Richards manufactured solution

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf, @sprintf
using Trixi

elixir = joinpath(dirname(@__DIR__), "elixirs",
                  "elixir_richards_manufactured_solution.jl")

# Convergence study options
time_steps = 4.0 ./ 2.0 .^ collect(0:3)
plot_ylims = (1.0e-12, 1.0e-3)

schemes = ((; name = "pressure_head", form = PressureHeadForm(), label = "Pressure head",
            color = 2, marker = :rect, markersize = 11.0),
           (; name = "mixed", form = MixedForm(), label = "Mixed", color = 1,
            marker = :circle, markersize = 7.0))

function solve_time_step(scheme, time_step)
    redirect_stdout(devnull) do
        Trixi.trixi_include(@__MODULE__, elixir; tspan = (0.0, 120.0),
                            initial_refinement_level = 8, polydeg = 3,
                            form = scheme.form, dt = time_step, adaptive = false,
                            reltol = 1.0e-9, abstol = 1.0e-11, saveat = Float64[],
                            run_simulation = true)
    end

    local_sol = Base.invokelatest(getproperty, @__MODULE__, :sol)
    local_analysis_callback = Base.invokelatest(getproperty, @__MODULE__,
                                                :analysis_callback)
    l2_error, linf_error = Base.invokelatest(local_analysis_callback, local_sol)
    return (; final_time = local_sol.t[end], nsteps = local_sol.destats.naccept,
            l2_error = only(l2_error), linf_error = only(linf_error),
            retcode = local_sol.retcode)
end

results = Dict{String, Any}()

for scheme in schemes
    final_times = Float64[]
    nsteps = Int[]
    l2_errors = Float64[]
    linf_errors = Float64[]
    retcodes = String[]

    for time_step in time_steps
        println("Running Richards manufactured solution $(scheme.name) with dt = $(time_step)")
        run = solve_time_step(scheme, time_step)
        push!(final_times, run.final_time)
        push!(nsteps, run.nsteps)
        push!(l2_errors, run.l2_error)
        push!(linf_errors, run.linf_error)
        push!(retcodes, string(run.retcode))
        println("retcode = $(run.retcode), steps = $(run.nsteps), L2 = $(run.l2_error), " *
                "Linf = $(run.linf_error)")
    end

    results[scheme.name] = (; time_steps = collect(time_steps), final_times, nsteps,
                            l2_errors, linf_errors, retcodes)
end

# Save one reproducible table for both formulations
plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
result_prefix = "richards_manufactured_solution_time_convergence"
analysis_path = joinpath(plots_dir, "$(result_prefix).dat")
open(analysis_path, "w") do io
    println(io,
            "form dt nsteps final_time l2_pressure_head_error " *
            "linf_pressure_head_error l2_order linf_order retcode")
    for scheme in schemes
        data = results[scheme.name]
        l2_orders = compute_eoc(data.l2_errors)
        linf_orders = compute_eoc(data.linf_errors)
        for i in eachindex(data.time_steps)
            @printf(io, "%s %.17e %d %.17e %.17e %.17e %.17e %.17e %s\n", scheme.name,
                    data.time_steps[i], data.nsteps[i], data.final_times[i],
                    data.l2_errors[i], data.linf_errors[i], l2_orders[i],
                    linf_orders[i], data.retcodes[i])
        end
    end
end
println("Saved Richards convergence table to: $(analysis_path)")

# Plot both formulations on one axis
HydroTrixi.set_serif_tex_theme!()
plot_font = HydroTrixi.DEFAULT_PLOT_FONT
time_step_ticks = Float64.(sort(unique(time_steps)))
time_step_labels = [@sprintf("%.4g", tick) for tick in time_step_ticks]
figure = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
axis = Axis(figure[1, 1]; xlabel = L"$\Delta t$ (s)",
            ylabel = LaTeXString("Pressure head error"), xlabelfont = plot_font,
            ylabelfont = plot_font, xticklabelfont = plot_font,
            yticklabelfont = plot_font, xscale = log10, yscale = log10,
            xticks = (time_step_ticks, time_step_labels))
axis.xminorgridvisible = false
axis.xminorticksvisible = false
isnothing(plot_ylims) || ylims!(axis, plot_ylims)

colors = Makie.wong_colors()
for scheme in schemes
    data = results[scheme.name]
    order = sortperm(data.time_steps)
    scatterlines!(axis, data.time_steps[order], data.l2_errors[order];
                  label = "$(scheme.label) L²", color = colors[scheme.color],
                  linestyle = :solid, linewidth = 1.8, marker = scheme.marker,
                  markersize = scheme.markersize)
    scatterlines!(axis, data.time_steps[order], data.linf_errors[order];
                  label = "$(scheme.label) L∞", color = colors[scheme.color],
                  linestyle = :dash, linewidth = 1.8, marker = scheme.marker,
                  markersize = scheme.markersize)
end
axislegend(axis; position = (:right, :bottom), labelsize = 14, font = plot_font)

output_path = joinpath(plots_dir, "$(result_prefix).pdf")
save(output_path, figure; px_per_unit = 1)
println("Saved Richards convergence plot to: $(output_path)")
