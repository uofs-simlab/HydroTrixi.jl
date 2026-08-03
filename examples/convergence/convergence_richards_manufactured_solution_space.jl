# Run and plot a fixed-mesh refinement study for the Richards manufactured solution

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__DIR__), "elixirs",
                  "elixir_richards_manufactured_solution.jl")

# Convergence study options
base_initial_refinement_level = 4
polydeg = 3

schemes = ((; name = "pressure_head", form = PressureHeadForm(), label = "Pressure head",
            color = 2, marker = :rect, markersize = 11.0),
           (; name = "mixed", form = MixedForm(), label = "Mixed", color = 1,
            marker = :circle, markersize = 7.0))

function solve_level(scheme, level)
    Trixi.trixi_include(@__MODULE__, elixir; tspan = (0.0, 120.0),
                        initial_refinement_level = level, polydeg = polydeg,
                        form = scheme.form, reltol = 1.0e-9, abstol = 1.0e-11,
                        dt = 1.0e-2, adaptive = true, saveat = Float64[],
                        run_simulation = true)

    local_sol = Base.invokelatest(getproperty, @__MODULE__, :sol)
    local_analysis_callback = Base.invokelatest(getproperty, @__MODULE__,
                                                :analysis_callback)
    l2_error, linf_error = Base.invokelatest(local_analysis_callback, local_sol)
    return (; final_time = local_sol.t[end], l2_error = only(l2_error),
            linf_error = only(linf_error), retcode = local_sol.retcode)
end

levels = collect(base_initial_refinement_level:(base_initial_refinement_level + 3))
results = Dict{String, Any}()

for scheme in schemes
    ndofs = (polydeg + 1) .* 2 .^ levels
    final_times = Float64[]
    l2_errors = Float64[]
    linf_errors = Float64[]

    for (level, ndof) in zip(levels, ndofs)
        println("Running Richards manufactured solution $(scheme.name) at level $(level)")
        run = solve_level(scheme, level)
        push!(final_times, run.final_time)
        push!(l2_errors, run.l2_error)
        push!(linf_errors, run.linf_error)
        println("retcode = $(run.retcode), ndofs = $(ndof), L2 = $(run.l2_error), " *
                "Linf = $(run.linf_error)")
    end

    results[scheme.name] = (; levels, ndofs, final_times, l2_errors, linf_errors)
end

# Save one reproducible table for both formulations
plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
result_prefix = "richards_manufactured_solution_space_convergence"
analysis_path = joinpath(plots_dir, "$(result_prefix).dat")
open(analysis_path, "w") do io
    println(io,
            "form level ndofs final_time l2_pressure_head_error " *
            "linf_pressure_head_error l2_order linf_order")
    for scheme in schemes
        data = results[scheme.name]
        l2_orders = compute_eoc(data.l2_errors)
        linf_orders = compute_eoc(data.linf_errors)
        for i in eachindex(data.levels)
            @printf(io, "%s %d %d %.17e %.17e %.17e %.17e %.17e\n", scheme.name,
                    data.levels[i], data.ndofs[i], data.final_times[i],
                    data.l2_errors[i], data.linf_errors[i], l2_orders[i],
                    linf_orders[i])
        end
    end
end
println("Saved Richards convergence table to: $(analysis_path)")

# Plot both formulations on one axis
HydroTrixi.set_serif_tex_theme!()
plot_font = HydroTrixi.DEFAULT_PLOT_FONT
all_ndofs = reduce(vcat, (results[scheme.name].ndofs for scheme in schemes))
figure = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
axis = Axis(figure[1, 1]; xlabel = LaTeXString("Degrees of freedom"),
            ylabel = LaTeXString("Pressure head error"), xlabelfont = plot_font,
            ylabelfont = plot_font, xticklabelfont = plot_font,
            yticklabelfont = plot_font, xscale = log10, yscale = log10,
            xticks = HydroTrixi.doubling_dof_ticks(all_ndofs; base = minimum(all_ndofs)))
axis.xminorgridvisible = false
axis.xminorticksvisible = false

colors = Makie.wong_colors()
for scheme in schemes
    data = results[scheme.name]
    scatterlines!(axis, data.ndofs, data.l2_errors; label = "$(scheme.label) L²",
                  color = colors[scheme.color], linestyle = :solid, linewidth = 1.8,
                  marker = scheme.marker, markersize = scheme.markersize)
    scatterlines!(axis, data.ndofs, data.linf_errors; label = "$(scheme.label) L∞",
                  color = colors[scheme.color], linestyle = :dash, linewidth = 1.8,
                  marker = scheme.marker, markersize = scheme.markersize)
end

reference_data = results[first(schemes).name]
y_ref = minimum(min(results[scheme.name].l2_errors[end - 1],
                    results[scheme.name].linf_errors[end - 1]) for scheme in schemes)
HydroTrixi.plot_bottom_triangle!(axis, reference_data.ndofs[end - 1],
                                 reference_data.ndofs[end], y_ref, polydeg + 1;
                                 trianglefontsize = 15, font = plot_font)
axislegend(axis; position = (:left, :bottom), labelsize = 14, font = plot_font)

output_path = joinpath(plots_dir, "$(result_prefix).pdf")
save(output_path, figure; px_per_unit = 1)
println("Saved Richards convergence plot to: $(output_path)")
