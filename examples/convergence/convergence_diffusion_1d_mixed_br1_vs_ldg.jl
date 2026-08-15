# Compare LDG and BR1 for the mixed Dirichlet-Neumann diffusion problem

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__DIR__), "elixirs",
                  "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")

# Convergence study options
diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 4

schemes = ((; name = "LDG", solver_parabolic = ParabolicFormulationLocalDG(),
            labels = (L"\mathrm{LDG}\ L^2", L"\mathrm{LDG}\ L^\infty"), color = 1),
           (; name = "BR1", solver_parabolic = ParabolicFormulationBassiRebay1(),
            labels = (L"\mathrm{BR1}\ L^2", L"\mathrm{BR1}\ L^\infty"), color = 2))

function run_scheme(scheme)
    _, errors = convergence_test(@__MODULE__, elixir, iterations;
                                 diffusivity = diffusivity,
                                 forcing_amplitude = forcing_amplitude,
                                 angular_frequency = 2 * pi * forcing_frequency,
                                 dirichlet_mean = 1.0, polydeg = polydeg,
                                 tspan = tspan,
                                 initial_refinement_level = base_initial_refinement_level,
                                 solver_parabolic = scheme.solver_parabolic)

    levels = collect(base_initial_refinement_level:
                     (base_initial_refinement_level + iterations - 1))
    return (; levels, ndofs = (polydeg + 1) .* 2 .^ levels,
            l2_errors = vec(errors[:l2]), linf_errors = vec(errors[:linf]))
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion convergence study with mixed Dirichlet-Neumann boundaries")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()

for scheme in schemes
    data = results[scheme.name]
    l2_eoc = compute_eoc(data.l2_errors)
    linf_eoc = compute_eoc(data.linf_errors)
    println("Scheme: $(scheme.name)")
    @printf("%-7s %-10s %-14s %-8s %-14s %-8s\n", "level", "ndofs", "L2 error",
            "EOC", "Linf error", "EOC")
    for i in eachindex(data.levels)
        @printf("%-7d %-10d %-14.6e %-8.3f %-14.6e %-8.3f\n", data.levels[i],
                data.ndofs[i], data.l2_errors[i], l2_eoc[i], data.linf_errors[i],
                linf_eoc[i])
    end
    println()
end

ldg = results["LDG"]
br1 = results["BR1"]
println("BR1/LDG error ratio by level")
@printf("%-7s %-12s %-12s\n", "level", "L2 ratio", "Linf ratio")
for i in eachindex(ldg.levels)
    @printf("%-7d %-12.4f %-12.4f\n", ldg.levels[i],
            br1.l2_errors[i]/ldg.l2_errors[i],
            br1.linf_errors[i]/ldg.linf_errors[i])
end

# Plot both parabolic formulations on one axis
series = map(schemes) do scheme
    data = results[scheme.name]
    (; x = data.ndofs, errors = (data.l2_errors, data.linf_errors),
     labels = scheme.labels, color = scheme.color)
end

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
output_path = joinpath(plots_dir, "diffusion_1d_mixed_br1_vs_ldg.pdf")
plot_convergence_1d(series; output_path, triangle_order = polydeg + 1,
                    trianglefontsize = 15)
println("Saved convergence plot to: $(output_path)")
