# Example: compare LDG and BR1 on the mixed Dirichlet-Neumann diffusion
# problem by running the same refinement study for both parabolic schemes.

using HydroTrixi
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__FILE__), "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")

diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
dirichlet_mean = 1.0
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 4
dt_factor = 0.01

schemes = ((; name = "LDG", solver_parabolic = ParabolicFormulationLocalDG()),
           (; name = "BR1", solver_parabolic = ParabolicFormulationBassiRebay1()))

function run_scheme(scheme)
    initial_refinement_level = base_initial_refinement_level
    _, errors = convergence_test(@__MODULE__,
                                 elixir,
                                 iterations;
                                 diffusivity = diffusivity,
                                 forcing_amplitude = forcing_amplitude,
                                 forcing_frequency = forcing_frequency,
                                 dirichlet_mean = dirichlet_mean,
                                 polydeg = polydeg,
                                 tspan = tspan,
                                 initial_refinement_level = initial_refinement_level,
                                 dt_factor = dt_factor,
                                 solver_parabolic = scheme.solver_parabolic)

    levels = initial_refinement_level:(initial_refinement_level + iterations - 1)
    ndofs = (polydeg + 1) .* (2 .^ collect(levels))

    return (; levels = collect(levels),
            ndofs,
            l2_errors = vec(errors[:l2]),
            linf_errors = vec(errors[:linf]))
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion equation convergence study (mixed Dirichlet-Neumann)")
println("Comparison: BR1 vs LDG parabolic schemes")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()

for scheme in schemes
    data = results[scheme.name]
    println("Scheme: $(scheme.name)")
    l2_eoc = compute_eoc(data.l2_errors)
    linf_eoc = compute_eoc(data.linf_errors)
    @printf("%-7s %-10s %-14s %-8s %-14s %-8s\n",
            "level",
            "ndofs",
            "L2 error",
            "EOC",
            "Linf error",
            "EOC")
    for i in eachindex(data.levels)
        @printf("%-7d %-10d %-14.6e ",
                data.levels[i],
                data.ndofs[i],
                data.l2_errors[i])
        isnan(l2_eoc[i]) ? @printf("%-8s ", "-") : @printf("%-8.3f ", l2_eoc[i])
        @printf("%-14.6e ", data.linf_errors[i])
        isnan(linf_eoc[i]) ? @printf("%-8s\n", "-") : @printf("%-8.3f\n", linf_eoc[i])
    end
    println()
end

ldg = results["LDG"]
br1 = results["BR1"]

println("BR1/LDG error ratio by level")
@printf("%-7s %-12s %-12s\n", "level", "L2 ratio", "Linf ratio")
for i in eachindex(ldg.levels)
    @printf("%-7d %-12.4f %-12.4f\n",
            ldg.levels[i],
            br1.l2_errors[i]/ldg.l2_errors[i],
            br1.linf_errors[i]/ldg.linf_errors[i])
end
