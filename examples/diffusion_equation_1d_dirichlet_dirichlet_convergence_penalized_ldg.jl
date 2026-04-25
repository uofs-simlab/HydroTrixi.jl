using DiffuSEM
using Printf: @printf
using SciMLBase
using Trixi

diffusivity = 0.5
polydeg = 3
tspan = (0.0, 0.25)
levels = collect(2:6)
dt_factor = 0.01
penalty_prefactor = 1.0

exact_solution(x, t) = exp(-diffusivity * pi^2 * t) * sinpi(x[1])
initial_condition(x, t, equations) = SVector(exact_solution(x, t))

function boundary_conditions_dirichlet(h)
    zero_dirichlet(x, t, equations) = SVector(0.0)
    return (;
            x_neg = BoundaryConditionDirichlet(zero_dirichlet),
            x_pos = BoundaryConditionDirichlet(zero_dirichlet),)
end

function boundary_conditions_penalized_dirichlet(h)
    zero_dirichlet(x, t, equations) = SVector(0.0)
    tau = penalty_prefactor * diffusivity * (polydeg + 1)^2 / h

    return (;
            x_neg = BoundaryConditionDirichletPenalty(zero_dirichlet; penalty = tau),
            x_pos = BoundaryConditionDirichletPenalty(zero_dirichlet; penalty = tau),)
end

schemes = ((; name = "LDG", boundary_builder = boundary_conditions_dirichlet),
           (;
            name = "LDG + boundary penalty",
            boundary_builder = boundary_conditions_penalized_dirichlet,))

function run_level(level, scheme)
    h = 1.0 / (2.0^level)
    mesh = TreeMesh((0.0,),
                    (1.0,),
                    initial_refinement_level = level,
                    n_cells_max = 30_000,
                    periodicity = false)
    equations = Trixi.LinearDiffusionEquation1D(diffusivity)
    solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)
    boundary_conditions = scheme.boundary_builder(h)

    semi = SemidiscretizationParabolic(mesh,
                                       equations,
                                       initial_condition,
                                       solver;
                                       boundary_conditions = boundary_conditions,
                                       solver_parabolic = ParabolicFormulationLocalDG(),)
    ode = semidiscretize(semi, tspan)

    dt = dt_factor * h^2
    sol = SciMLBase.solve(ode,
                          default_algorithm();
                          dt = dt,
                          adaptive = false,
                          save_everystep = false,)

    errors = AnalysisCallback(semi, interval = typemax(Int))(sol)
    return length(sol.u[end]), errors.l2[1], errors.linf[1], dt
end

function run_scheme(scheme)
    ndofs = Int[]
    l2_errors = Float64[]
    linf_errors = Float64[]
    dts = Float64[]

    for level in levels
        dofs, l2, linf, dt = run_level(level, scheme)
        push!(ndofs, dofs)
        push!(l2_errors, l2)
        push!(linf_errors, linf)
        push!(dts, dt)
    end

    return (;
            ndofs,
            l2_errors,
            linf_errors,
            dts,
            l2_eoc = compute_eoc(l2_errors),
            linf_eoc = compute_eoc(linf_errors),)
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion equation convergence study (double Dirichlet)")
println("u(x,t) = exp(-nu*pi^2*t) * sin(pi*x), nu = $(diffusivity)")
println("polydeg = $(polydeg), tspan = $(tspan), penalty_prefactor = $(penalty_prefactor)")
println()

for scheme in schemes
    data = results[scheme.name]
    println("Scheme: $(scheme.name)")
    @printf("%-7s %-10s %-11s %-14s %-8s %-14s %-8s\n",
            "level",
            "ndofs",
            "dt",
            "L2 error",
            "EOC",
            "Linf error",
            "EOC")
    for i in eachindex(levels)
        @printf("%-7d %-10d %-11.3e %-14.6e ",
                levels[i],
                data.ndofs[i],
                data.dts[i],
                data.l2_errors[i])
        isnan(data.l2_eoc[i]) ? @printf("%-8s ", "-") : @printf("%-8.3f ", data.l2_eoc[i])
        @printf("%-14.6e ", data.linf_errors[i])
        isnan(data.linf_eoc[i]) ? @printf("%-8s\n", "-") :
        @printf("%-8.3f\n", data.linf_eoc[i])
    end
    println()
end
