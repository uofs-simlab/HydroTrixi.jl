# Render an AMR animation of the Celia (1990) pressure-head form
# The animation routine advances the ODE problem so each frame uses its mesh

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

final_time = 360.0
saveat = 0.0:2.0:final_time

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
result_prefix = "richards_celia_1990_pressure_head_form_amr_interval30_maxlevel6_t" *
                "$(round(Int, final_time))"
analysis_filename = "$(result_prefix)_analysis.dat"
analysis_path = joinpath(plots_dir, analysis_filename)

Trixi.trixi_include(joinpath(@__DIR__, "elixir_richards_celia_1990_amr.jl");
                    final_time = final_time,
                    form = PressureHeadForm(),
                    run_simulation = false,
                    saveat = saveat,
                    amr_interval = 30,
                    max_level = 6,
                    save_analysis = true,
                    output_directory = plots_dir,
                    analysis_filename = analysis_filename,
                    reltol = 1.0e-7,
                    abstol = 1.0e-11)

if !@isdefined(animation_format)
    animation_format = "mp4"
end
if !@isdefined(plot_format)
    plot_format = "pdf"
end

output_path = joinpath(plots_dir, "$(result_prefix).$(animation_format)")
mass_bias_output_path = joinpath(plots_dir, "$(result_prefix)_mass_bias.$(plot_format)")

animate_solution_1d(ode;
                    callback = callbacks,
                    dt = dt,
                    adaptive = adaptive,
                    reltol = reltol,
                    abstol = abstol,
                    saveat = saveat,
                    component = 1,
                    xlabel = L"$z$ (m)",
                    ylabel = L"$\psi$ (m)",
                    ylims = (-0.65, -0.15),
                    show_element_boundaries = true,
                    output_path = output_path,
                    framerate = 30)

plot_mass_bias(analysis_path;
               output_path = mass_bias_output_path,
               xlabel = L"$t$ (s)",
               ylabel = L"$\epsilon_b$ (m)",
               xticks = 0.0:60.0:final_time)

println("Saved Celia pressure head AMR animation to: $(output_path)")
println("Saved Celia pressure head AMR mass-bias plot to: $(mass_bias_output_path)")
println("Saved Celia pressure head AMR analysis data to: $(analysis_path)")
