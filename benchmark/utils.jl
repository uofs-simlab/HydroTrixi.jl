using BenchmarkTools
using Dates
using InteractiveUtils: versioninfo
using Statistics: median

# Print runtime-generated metadata needed to reproduce a benchmark
function print_benchmark_environment(packages; repository_directory)
    revision = readchomp(`git -C $repository_directory rev-parse HEAD`)
    working_tree_status = read(`git -C $repository_directory status --porcelain`, String)

    println("Date: $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))")
    println("Repository revision: $(revision)")
    println("Working tree dirty: $(!isempty(working_tree_status))")
    println()

    println("System information")
    versioninfo()
    println()

    println("Package versions")
    for package in packages
        println("  $(nameof(package)).jl = $(pkgversion(package))")
    end
    println()
    return nothing
end

# Collect independent evaluations after the automatic BenchmarkTools.jl warm-up
function run_benchmark(function_to_benchmark, arguments...; samples, seconds)
    target = () -> function_to_benchmark(arguments...)
    benchmark = @benchmarkable $target()
    trial = BenchmarkTools.run(benchmark; samples, evals = 1, seconds,
                               gcsample = true)
    length(trial) == samples ||
        error("Benchmark collected only $(length(trial)) of $(samples) samples")
    return trial
end

# Print a standard BenchmarkTools.jl trial summary
function print_benchmark_summary(label, trial)
    minimum_estimate = minimum(trial)
    median_estimate = median(trial)
    maximum_estimate = maximum(trial)
    minimum_time = BenchmarkTools.time(minimum_estimate) / 1.0e9
    median_time = BenchmarkTools.time(median_estimate) / 1.0e9
    maximum_time = BenchmarkTools.time(maximum_estimate) / 1.0e9
    median_gc_time = BenchmarkTools.gctime(median_estimate) / 1.0e9
    allocated_mebibytes = BenchmarkTools.memory(trial) / 2^20

    println("  $(label)")
    println("    median time = $(round(median_time; sigdigits = 5)) s")
    println("    time range = $(round(minimum_time; sigdigits = 5))-" *
            "$(round(maximum_time; sigdigits = 5)) s")
    println("    memory estimate = " *
            "$(round(allocated_mebibytes; sigdigits = 5)) MiB")
    println("    median GC time = $(round(median_gc_time; sigdigits = 5)) s")

    return (; minimum_time, median_time, maximum_time, median_gc_time,
            allocated_mebibytes)
end
