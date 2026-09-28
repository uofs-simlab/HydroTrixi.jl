# Run all six studies or a selected subset; see README.md.

module RichardsConvergence

using HydroTrixi, Trixi, SciMLBase
using Dates, LinearAlgebra

include("plot_richards_convergence.jl")

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const ELIXIR = joinpath(ROOT, "examples", "elixirs",
                        "elixir_richards_manufactured_solution.jl")
const FORMS = (("pressure_head", PressureHeadForm()), ("mixed", MixedForm()))
const STUDIES = [(; bc, kind, N) for bc in ("DD", "DN")
                 for (kind, N) in (("space", 3), ("space", 4), ("time", 3))]

# Keep the elixir's global assignments separate from the convergence driver.
module ManufacturedSolutionElixir end

function solve_case(bc, form, level, N, dt)
    final_time = 120.0
    problem_options = (;)
    if bc == "DN"
        base = HydrologicProblemRichardsManufacturedSolution()
        boundaries = (; x_neg = base.boundary_conditions.x_neg,
                        x_pos = BoundaryConditionNeumann(HydroTrixi.richards_manufactured_right_boundary_flux))
        problem = HydrologicProblem(; equations = base.equations,
            state_to_evolved = base.state_to_evolved, evolved_to_state = base.evolved_to_state,
            initial_condition = base.initial_condition, boundary_conditions = boundaries,
            source_terms = base.source_terms, domain = base.domain, tspan = base.tspan)
        problem_options = (; problem)
    end
    # Prescribe the grid to avoid a roundoff-sized extra step for decimal dt.
    nsteps = round(Int, final_time / dt)
    time_grid = range(0.0, final_time; length = nsteps + 1)[2:end]
    simulation = Trixi.trixi_include(ManufacturedSolutionElixir, ELIXIR;
                                     problem_options..., initial_refinement_level = level,
                                     polydeg = N, form, dt, adaptive = false,
                                     callback = nothing,
                                     solve_options = (; tstops = time_grid, dense = false))
    (; sol, analysis_callback) = simulation

    l2, linf = analysis_callback(sol)
    valid = SciMLBase.successful_retcode(sol)
    return (; l2 = only(l2), linf = only(linf), final_time = last(sol.t),
            steps = sol.stats.naccept, retcode = string(sol.retcode), valid)
end

function run_convergence(studies = STUDIES)
    if !all(study -> study in STUDIES, studies)
        error("Unknown convergence study")
    end
    BLAS.set_num_threads(1)

    # Give each run its own directory, preserving earlier and partial results.
    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output = joinpath(ROOT, "plots", "richards_convergence", id)
    mkpath(dirname(output))
    mkdir(output)
    mkdir(joinpath(output, "data"))

    for study in studies
        (; bc, kind, N) = study
        boundary_name = bc == "DD" ? "dirichlet_dirichlet" : "dirichlet_neumann"
        name = "richards_manufactured_solution_$(boundary_name)_$(kind)_N$(N)"
        println("Running $name")
        configs = kind == "space" ? [(level, 0.1) for level in 4:8] :
                                   [(12, dt) for dt in (4.0, 2.0, 1.0, 0.5, 0.25)]

        # Write each result immediately so partial runs remain available.
        open(joinpath(output, "data", "$name.dat"), "w") do data
            println(data, "form level N elements pressure_dofs total_dofs delta_z_m dt_s final_time_s accepted_steps l2_m linf_m l2_order linf_order retcode status")
            flush(data)
            for (form_name, form) in FORMS
                previous = nothing
                for (level, dt) in configs
                    result = redirect_stdout(devnull) do
                        solve_case(bc, form, level, N, dt)
                    end
                    p2 = isnothing(previous) ? NaN : log2(previous.l2 / result.l2)
                    pinf = isnothing(previous) ? NaN : log2(previous.linf / result.linf)
                    cells = 2^level
                    dofs = cells * (N + 1)
                    status = result.valid ? "success" : "failed"
                    println(data, join((form_name, level, N, cells, dofs,
                        dofs * (form isa MixedForm ? 2 : 1), 0.2 / (cells * N), dt,
                        result.final_time, result.steps, result.l2, result.linf, p2, pinf,
                        result.retcode, status), ' '))
                    flush(data)
                    if !result.valid
                        error("Failed $name $form_name level=$level dt=$dt s; see $output")
                    end
                    previous = result
                end
            end
        end
    end
    RichardsConvergencePlots.plot_convergence(output)
    println("Results: $output")
    return output
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) in (0, 3))
        error("Usage: julia --project=run $(@__FILE__) [DD|DN space|time N]")
    end
    studies = isempty(ARGS) ? RichardsConvergence.STUDIES :
              [(bc = ARGS[1], kind = ARGS[2], N = parse(Int, ARGS[3]))]
    RichardsConvergence.run_convergence(studies)
end
