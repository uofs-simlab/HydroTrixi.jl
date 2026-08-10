# Benchmark the fixed-mesh Celia problem in mixed and pressure-head forms

using HydroTrixi
using SciMLBase: solve, successful_retcode
using Trixi

include("utils.jl")

# Benchmark configuration
final_time = 360.0
time_step = 1.0e-2
samples = 5
seconds_per_trial = 120.0

schemes = ((; name = :pressure_head, label = "Pressure head",
            form = PressureHeadForm()),
           (; name = :mixed, label = "Mixed", form = MixedForm()))

# Run the complete problem setup and solve within the timing boundary
function run_case(form, case_final_time, time_step)
    problem = HydrologicProblemCelia1990(tspan = (0.0, case_final_time))
    mesh = TreeMesh(problem.domain..., initial_refinement_level = 5,
                    periodicity = false)
    solver = DGSEM(polydeg = 3)
    solver_parabolic = ParabolicFormulationLocalDG()
    passive_variables = PassiveVariablesBoundaryFlux1D()
    semidiscretization = SemidiscretizationImplicit(mesh, problem, solver;
                                                    solver_parabolic,
                                                    passive_variables,
                                                    form)
    ode = semidiscretize(semidiscretization, problem.tspan;
                         jacobian = SparseJacobian())

    solution = solve(ode, default_algorithm(semidiscretization);
                     dt = time_step, adaptive = false, saveat = Float64[],
                     ode_default_options()..., maxiters = typemax(Int))

    return (; successful = successful_retcode(solution),
            final_time = solution.t[end],
            steps = solution.destats.naccept,
            rejected_steps = solution.destats.nreject,
            rhs_calls = solution.destats.nf,
            spatial_degrees_of_freedom = Trixi.ndofs(semidiscretization),
            unknowns = length(ode.u0),
            elements = length(Trixi.leaf_cells(mesh.tree)),
            mass_bias = HydroTrixi.mass_bias(solution.u[end], semidiscretization),
            autodiff = string(typeof(solution.alg.autodiff)))
end

print_benchmark_environment((BenchmarkTools, HydroTrixi, Trixi);
                            repository_directory = dirname(@__DIR__))

println("Fixed-mesh Celia benchmark")
println("  tspan = (0.0, $(final_time))")
println("  dt = $(time_step), adaptive = false, AMR = false")
println("  samples = $(samples), timing boundary = setup + solve + diagnostics")

# Validate and benchmark each formulation
println("Results")
median_times = Dict{Symbol, Float64}()
for scheme in schemes
    reference_result = run_case(scheme.form, final_time, time_step)
    reference_result.successful || error("$(scheme.label) benchmark solve failed")
    reference_result.final_time == final_time ||
        error("$(scheme.label) benchmark solve stopped before its final time")
    occursin("AutoSparse", reference_result.autodiff) ||
        error("$(scheme.label) solve did not prepare sparse automatic differentiation")

    trial = run_benchmark(run_case, scheme.form, final_time, time_step;
                          samples, seconds = seconds_per_trial)
    summary = print_benchmark_summary(scheme.label, trial)
    median_times[scheme.name] = summary.median_time

    println("    steps = $(reference_result.steps), rejected = " *
            "$(reference_result.rejected_steps), RHS calls = " *
            "$(reference_result.rhs_calls)")
    println("    elements = $(reference_result.elements), spatial DOFs per field = " *
            "$(reference_result.spatial_degrees_of_freedom)")
    println("    ODE/DAE unknowns = $(reference_result.unknowns)")
    println("    absolute mass bias = " *
            "$(round(abs(reference_result.mass_bias); sigdigits = 5))")
    println("    autodiff = $(reference_result.autodiff)")
end

runtime_ratio = median_times[:pressure_head] / median_times[:mixed]
println("  pressure head to mixed median runtime ratio = " *
        "$(round(runtime_ratio; sigdigits = 5))")
