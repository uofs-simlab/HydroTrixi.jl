# Render an AMR animation of the closed-column Richards-equation redistribution profile
# The animation routine advances the ODE problem so each frame uses its mesh

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi: trixi_include

final_time = 86400.0

saveat = range(0.0, final_time; length = 181)

trixi_include(joinpath(@__DIR__, "elixir_richards_closed_column_amr.jl");
              final_time = final_time,
              run_simulation = false,
              saveat = saveat)

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(animation_format)
    animation_format = "mp4"
end

output_path = joinpath(plots_dir,
                       "richards_closed_column_amr_pressure_head_t$(round(Int, final_time)).$(animation_format)")

animate_solution_1d(ode;
                    callback = callbacks,
                    dt = dt,
                    adaptive = adaptive,
                    saveat = saveat,
                    component = 2,
                    xlabel = L"$z$ (m)",
                    ylabel = L"$\psi$ (m)",
                    ylims = (-1.05, -0.55),
                    show_element_boundaries = true,
                    output_path = output_path,
                    framerate = 30)

println("Saved closed-column AMR pressure head animation to: $(output_path)")
