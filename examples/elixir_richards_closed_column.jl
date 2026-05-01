# Example: Trixi-style elixir for a closed-column Richards-equation
# redistribution problem with zero normal flux at both ends. The initial
# pressure head profile satisfies dpsi/dz = 1 at the boundaries, so the
# physical flux K(psi) * (dpsi/dz - 1) vanishes initially as well. This makes
# the setup useful for assessing discrete mass conservation via the Richards
# water-content analysis quantities.

using HydroTrixi
using SciMLBase
using Trixi

if !@isdefined(final_time)
    final_time = 360.0
end

initial_refinement_level = 5
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()
adaptive = true
dt = 1.0e-2
if !@isdefined(saveat)
    saveat = 0.0:10.0:final_time
end
save_everystep = false

problem = HydrologicProblemRichardsClosedColumn(tspan = (0.0, final_time))
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
analysis_callback = AnalysisCallback(semi,
                                     interval = 10,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = (HydroTrixi.water_content,
                                                                 HydroTrixi.water_content_timederivative))
ode = semidiscretize(semi, tspan)

sol = solve(ode,
            default_algorithm(semi);
            callback = analysis_callback,
            dt = dt,
            adaptive = adaptive,
            saveat = saveat,
            save_everystep = save_everystep,
            maxiters = typemax(Int))
