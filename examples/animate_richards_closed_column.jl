# Render an animation of the closed-column Richards-equation redistribution
# profile. Re-runs the elixir with a denser `saveat` so the solution carries
# enough frames for a smooth animation, then animates the pressure head block
# of the constitutive ODE state (component = 2).

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi: trixi_include

if !@isdefined(final_time)
    final_time = 1440.0
end

if !@isdefined(saveat)
    saveat = 0.0:8.0:final_time
end

trixi_include(joinpath(@__DIR__, "elixir_richards_closed_column.jl");
              final_time = final_time,
              saveat = saveat)

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(animation_format)
    animation_format = "mp4"
end

output_path = joinpath(plots_dir,
                       "richards_closed_column_pressure_head_t$(round(Int, final_time)).$(animation_format)")

animate_solution_1d(sol;
                    component = 2,
                    xlabel = L"$z$ (m)",
                    ylabel = L"$\psi$ (m)",
                    ylims = (-1.05, -0.55),
                    output_path = output_path,
                    framerate = 30)

println("Saved closed-column ψ(z,t) animation to: $(output_path)")
