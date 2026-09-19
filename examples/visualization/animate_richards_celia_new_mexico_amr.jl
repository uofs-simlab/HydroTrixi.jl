# Animate the spatially adaptive Celia New Mexico benchmark

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

# Simulation and visualization options
form = MixedForm()
final_time = 86_400.0
amr_interval = parse(Int, get(ENV, "HYDROTRIXI_CELIA_NM_AMR_INTERVAL", "10"))
amr_base_level = parse(Int, get(ENV, "HYDROTRIXI_CELIA_NM_AMR_BASE_LEVEL", "2"))

pressure_head_form = form isa PressureHeadForm
form_name = pressure_head_form ? "pressure_head" : "mixed"
component = pressure_head_form ? 1 : 2
result_prefix = "richards_celia_new_mexico_$(form_name)_amr_interval$(amr_interval)_base$(amr_base_level)_t$(round(Int, final_time))"

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
analysis_filename = "$(result_prefix)_analysis.dat"
analysis_path = joinpath(plots_dir, analysis_filename)

Trixi.trixi_include(@__MODULE__,
                    joinpath(dirname(@__DIR__), "elixirs",
                             "elixir_richards_celia_new_mexico.jl");
                    tspan = (0.0, final_time), form = form, analysis_interval = 200,
                    amr = true, amr_interval = amr_interval,
                    amr_base_level = amr_base_level,
                    run_simulation = false, save_analysis = true,
                    output_directory = plots_dir, analysis_filename = analysis_filename)

animation_path = joinpath(plots_dir, "$(result_prefix).mp4")
mass_bias_path = joinpath(plots_dir, "$(result_prefix)_mass_bias.pdf")

animate_solution_1d(ode; callback = callbacks, dt = 1.0e-2, adaptive = true,
                    reltol = 1.0e-7, abstol = 1.0e-11,
                    saveat = range(0.0, final_time; length = 181),
                    component = component, xlabel = L"$z$ (m)", ylabel = L"$\psi$ (m)",
                    ylims = (-10.5, -0.5), show_element_boundaries = true,
                    output_path = animation_path, framerate = 30)

plot_mass_bias(analysis_path; output_path = mass_bias_path, xlabel = L"$t$ (s)",
               ylabel = L"$\epsilon_b$ (m)",
               xticks = range(0.0, final_time; length = 7))

println("Saved Celia New Mexico $(form_name) AMR animation to: $(animation_path)")
println("Saved Celia New Mexico $(form_name) AMR mass-bias plot to: $(mass_bias_path)")
println("Saved Celia New Mexico $(form_name) AMR analysis data to: $(analysis_path)")
