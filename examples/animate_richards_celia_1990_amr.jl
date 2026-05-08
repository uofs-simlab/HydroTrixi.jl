# Render an AMR animation of the Celia (1990) Richards-equation infiltration
# profile. The animation routine advances the ODE problem so each frame is
# plotted against the mesh used at that time.

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi: trixi_include

saveat = 0.0:2.0:360.0

trixi_include(joinpath(@__DIR__, "elixir_richards_celia_1990_amr.jl");
              run_simulation = false,
              saveat = saveat)

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(animation_format)
    animation_format = "mp4"
end

output_path = joinpath(plots_dir,
                       "richards_celia_1990_amr_pressure_head.$(animation_format)")

animate_solution_1d(ode;
                    callback = callbacks,
                    dt = dt,
                    adaptive = adaptive,
                    saveat = saveat,
                    component = 2,
                    xlabel = L"$z$ (m)",
                    ylabel = L"$\psi$ (m)",
                    ylims = (-0.65, -0.15),
                    show_element_boundaries = true,
                    output_path = output_path,
                    framerate = 30)

println("Saved Celia AMR pressure-head animation to: $(output_path)")
@show last(saveat)
