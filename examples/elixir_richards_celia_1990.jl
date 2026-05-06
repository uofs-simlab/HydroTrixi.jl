# Example: Trixi-style elixir for the Celia et al. (1990) Richards equation infiltration
# problem. The problem definition lives in `HydrologicProblemCelia1990()`
# (Haverkamp soil model, Dirichlet at both ends, depth z in m on [0, 0.4],
# t in seconds on [0, 360]).

using HydroTrixi
using SciMLBase
using Trixi

initial_refinement_level = 5    # 32 elements over 0.4 m → Δz ≈ 0.0125 m
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()
adaptive = true
dt = 1.0e-2                     # initial step; Rodas5P adapts from here
saveat = Float64[]
save_everystep = false

problem = HydrologicProblemCelia1990()
(; domain, tspan) = problem
lower, upper = domain

mesh = TreeMesh(lower,
                upper,
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max,
                periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = solver_parabolic)
ode = semidiscretize(semi, tspan)

sol = solve(ode,
            default_algorithm(semi);
            dt = dt,
            adaptive = adaptive,
            saveat = saveat,
            save_everystep = save_everystep,
            maxiters = typemax(Int))
