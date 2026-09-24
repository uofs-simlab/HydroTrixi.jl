# Run all eight conservation studies or one selected study; see README.md.

module RichardsConservation

using HydroTrixi, Trixi, SciMLBase
using Dates, LinearAlgebra

include("plot_richards_conservation.jl")

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCHMARKS = (
    (name = "haverkamp", prefix = "richards_celia_haverkamp", final_time = 360.0,
     elixir = joinpath(ROOT, "examples", "elixirs",
                      "elixir_richards_celia_haverkamp.jl")),
    (name = "new_mexico", prefix = "richards_celia_new_mexico",
     final_time = 86_400.0,
     elixir = joinpath(ROOT, "examples", "elixirs",
                      "elixir_richards_celia_new_mexico.jl")),
)
const CASES = (
    (name = "mixed", suffix = "mixed", form = MixedForm()),
    (name = "mixed_pressure_head_transfer", suffix = "mixed_pressure_head_transfer",
     form = MixedForm(; transfer_state = true)),
    (name = "pressure_head", suffix = "pressure_head_pressure_head_transfer",
     form = PressureHeadForm()),
    (name = "pressure_head_water_content_transfer",
     suffix = "pressure_head_water_content_transfer",
     form = PressureHeadForm(; transfer_variables = water_content)),
)
const RELTOLS = RichardsConservationPlots.RELTOLS
const STUDIES = [(benchmark = benchmark.name, case = case.name)
                 for benchmark in BENCHMARKS for case in CASES]

# Keep the elixirs' global assignments separate from the conservation driver. Predeclare
# the bindings read by compiled driver functions before `trixi_include` assigns them.
module InfiltrationElixir
ode = nothing
amr_callback = nothing
analysis_callback = nothing
end

function select_config(configs, name)
    index = findfirst(config -> config.name == name, configs)
    if isnothing(index)
        error("Unknown conservation-study option: $name")
    end
    return configs[index]
end

function study_name(study)
    benchmark = select_config(BENCHMARKS, study.benchmark)
    case = select_config(CASES, study.case)
    return "$(benchmark.prefix)_$(case.suffix)"
end

function read_analysis(path)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    if isempty(lines) || !startswith(strip(first(lines)), "#")
        error("Invalid analysis table: $path")
    end
    columns = split(strip(first(lines))[2:end])
    indices = Dict(column => index for (index, column) in pairs(columns))
    required = ("timestep", "time", "dt", "mass_balance")
    if !all(haskey(indices, column) for column in required)
        error("Missing columns in analysis table: $path")
    end

    rows = map(lines[2:end]) do line
        values = split(line)
        (; step = parse(Int, values[indices["timestep"]]),
         time = parse(Float64, values[indices["time"]]),
         dt = parse(Float64, values[indices["dt"]]),
         balance = parse(Float64, values[indices["mass_balance"]]))
    end
    if isempty(rows)
        error("Analysis table contains no samples: $path")
    end
    initial_balance = first(rows).balance
    return [(; row.step, row.time, row.dt, bias = row.balance - initial_balance)
            for row in rows]
end

function validate_result(sol, rows, final_time, partial_path)
    accepted_steps = sol.stats.naccept
    valid = SciMLBase.successful_retcode(sol) && last(sol.t) == final_time &&
            length(rows) == accepted_steps + 1 &&
            getproperty.(rows, :step) == collect(0:accepted_steps) &&
            last(rows).time == final_time &&
            all(row -> all(isfinite, (row.time, row.dt, row.bias)), rows)
    if !valid
        error("Invalid conservation result; partial data retained at $partial_path")
    end
    return nothing
end

function write_mass_bias_table(path, rows)
    open(path, "w") do io
        println(io, "#accepted_step time_s dt_s mass_bias_m")
        for row in rows
            println(io, join((row.step, row.time, row.dt, row.bias), ' '))
        end
    end
    return path
end

function solve_case(study, tolerance, data_directory)
    benchmark = select_config(BENCHMARKS, study.benchmark)
    case = select_config(CASES, study.case)
    name = study_name(study)
    partial_path = joinpath(data_directory,
                            ".$(name)_rtol$(tolerance.tag).partial")
    final_path = joinpath(data_directory,
                          RichardsConservationPlots.data_filename(name, tolerance))

    Trixi.trixi_include(InfiltrationElixir, benchmark.elixir;
                        tspan = (0.0, benchmark.final_time), polydeg = 3,
                        initial_refinement_level = 6, form = case.form,
                        amr = true, amr_interval = 10, base_level = 2,
                        coarsen_threshold = 0.003, max_level = 10,
                        refine_threshold = 0.03, adapt_initial_condition = true,
                        adapt_initial_condition_only_refine = true,
                        analysis_interval = 1, save_analysis = true,
                        output_directory = data_directory,
                        analysis_filename = basename(partial_path),
                        run_simulation = false)

    ode = InfiltrationElixir.ode
    callbacks = CallbackSet(InfiltrationElixir.amr_callback,
                            InfiltrationElixir.analysis_callback)
    pressure_head_form = case.form isa PressureHeadForm
    error_control_block, error_control_mapping = if pressure_head_form
        (state_variable_block, water_content)
    else
        (evolved_variable_block, nothing)
    end
    sol = redirect_stdout(devnull) do
        solve_implicit(ode; dt = 1.0e-2, dtmin = 0.0, adaptive = true,
                       reltol = tolerance.value, abstol = 1.0e-11,
                       saveat = Float64[], dense = false,
                       error_control_block, error_control_mapping,
                       isoutofdomain = pressure_head_out_of_domain,
                       callback = callbacks)
    end

    rows = read_analysis(partial_path)
    validate_result(sol, rows, benchmark.final_time, partial_path)
    write_mass_bias_table(final_path, rows)
    rm(partial_path)
    return final_path
end

function run_study(study, output)
    for tolerance in RELTOLS
        println("Running $(study_name(study)) with reltol=$(tolerance.tag)")
        solve_case(study, tolerance, joinpath(output, "data"))
    end
    return nothing
end

function run_conservation(studies = STUDIES)
    if isempty(studies) || !all(study -> study in STUDIES, studies) ||
       !allunique(studies)
        error("Invalid study selection")
    end
    BLAS.set_num_threads(1)
    if Threads.nthreads() != 1
        error("Run with JULIA_NUM_THREADS=1")
    end

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output = joinpath(ROOT, "plots", "richards_conservation", id)
    mkpath(dirname(output))
    mkdir(output)
    mkdir(joinpath(output, "data"))

    for study in studies
        run_study(study, output)
    end
    RichardsConservationPlots.plot_conservation(output)
    println("Results: $output")
    return output
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) in (0, 2))
        error("Usage: julia --project=run $(@__FILE__) " *
              "[haverkamp|new_mexico " *
              "mixed|mixed_pressure_head_transfer|pressure_head|" *
              "pressure_head_water_content_transfer]")
    end
    studies = if isempty(ARGS)
        RichardsConservation.STUDIES
    else
        [(benchmark = ARGS[1], case = ARGS[2])]
    end
    RichardsConservation.run_conservation(studies)
end
