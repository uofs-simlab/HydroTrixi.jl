using HydroTrixi
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__FILE__), "elixir_diffusion_1d_dirichlet_dirichlet.jl")

diffusivity = 0.5
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 5
dt_factor = 0.01
penalty_prefactor = 1.0

schemes = ((; name = "LDG", use_boundary_penalty = false),
           (; name = "LDG + boundary penalty", use_boundary_penalty = true))

function run_scheme(scheme)
    initial_refinement_level = base_initial_refinement_level
    _, errors = convergence_test(@__MODULE__, elixir, iterations; diffusivity = diffusivity,
                                 polydeg = polydeg, tspan = tspan,
                                 initial_refinement_level = initial_refinement_level,
                                 dt_factor = dt_factor,
                                 solver_parabolic = ParabolicFormulationLocalDG(),
                                 use_boundary_penalty = scheme.use_boundary_penalty,
                                 penalty_prefactor = penalty_prefactor)

    levels = initial_refinement_level:(initial_refinement_level + iterations - 1)

    return (; levels = collect(levels), ndofs = (polydeg + 1) .* (2 .^ collect(levels)),
            dts = dt_factor .* (1.0 ./ (2.0 .^ collect(levels))) .^ 2,
            l2_errors = vec(errors[:l2]), linf_errors = vec(errors[:linf]))
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion equation convergence study (double Dirichlet)")
println("u(x,t) = exp(-nu*pi^2*t) * sin(pi*x), nu = $(diffusivity)")
println("polydeg = $(polydeg), tspan = $(tspan), penalty_prefactor = $(penalty_prefactor)")
println()

for scheme in schemes
    data = results[scheme.name]
    println("Scheme: $(scheme.name)")
    l2_eoc = compute_eoc(data.l2_errors)
    linf_eoc = compute_eoc(data.linf_errors)
    @printf("%-7s %-10s %-11s %-14s %-8s %-14s %-8s\n", "level", "ndofs", "dt", "L2 error",
            "EOC", "Linf error", "EOC")
    for i in eachindex(data.levels)
        @printf("%-7d %-10d %-11.3e %-14.6e ", data.levels[i], data.ndofs[i], data.dts[i],
                data.l2_errors[i])
        isnan(l2_eoc[i]) ? @printf("%-8s ", "-") : @printf("%-8.3f ", l2_eoc[i])
        @printf("%-14.6e ", data.linf_errors[i])
        isnan(linf_eoc[i]) ? @printf("%-8s\n", "-") :
        @printf("%-8.3f\n", linf_eoc[i])
    end
    println()
end
