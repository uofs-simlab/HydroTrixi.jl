# Animate the closed-column Richards problem with a fixed or adaptive mesh

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

# Simulation and visualization options
amr = false
final_time = amr ? 86400.0 : 1440.0
saveat = range(0.0, final_time; length = 181)

Trixi.trixi_include(@__MODULE__,
                    joinpath(dirname(@__DIR__), "elixirs",
                             "elixir_richards_closed_column.jl");
                    tspan = (0.0, final_time), amr = amr, run_simulation = !amr,
                    saveat = saveat)

mesh_name = amr ? "amr" : "fixed"
plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
output_path = joinpath(plots_dir,
                       "richards_closed_column_$(mesh_name)_t$(round(Int, final_time)).mp4")

animation_options = (; component = 2, xlabel = L"$z$ (m)", ylabel = L"$\psi$ (m)",
                     ylims = (-1.05, -0.55), output_path = output_path, framerate = 30)
if amr
    animate_solution_1d(ode; callback = callbacks, dt = 1.0e-2, adaptive = true,
                        saveat = saveat, show_element_boundaries = true,
                        animation_options...)
else
    animate_solution_1d(sol; animation_options...)
end

println("Saved closed-column $(mesh_name)-mesh animation to: $(output_path)")
