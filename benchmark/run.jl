using Dates
using Pkg

# Prepare the benchmark environment against the local HydroTrixi.jl checkout
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
Pkg.instantiate()

benchmark_filename = isempty(ARGS) ? "benchmark_richards_celia_1990_fixed.jl" : only(ARGS)
benchmark_path = isabspath(benchmark_filename) ? benchmark_filename :
                 joinpath(@__DIR__, benchmark_filename)
result_filename = splitext(basename(benchmark_path))[1] * "_" *
                  Dates.format(today(), dateformat"yyyy-mm-dd") * ".txt"
result_path = joinpath(@__DIR__, "results", result_filename)

# Archive exactly the benchmark standard output and echo it after completion
mkpath(dirname(result_path))
open(result_path, "w") do result_stream
    redirect_stdout(result_stream) do
        include(benchmark_path)
    end
end
print(read(result_path, String))
