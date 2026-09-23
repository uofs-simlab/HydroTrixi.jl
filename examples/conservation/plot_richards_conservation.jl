# Replot one or more saved conservation-study directories; never runs a simulation.

module RichardsConservationPlots

using CairoMakie, HydroTrixi, LaTeXStrings
using Dates, Printf

const RTOL_FONT = assetpath("fonts", "DejaVuSansMono.ttf")
const RELTOLS = ((value = 1.0e-5, tag = "1e-5",
                  label = rich(rich("rtol", font = RTOL_FONT), " = 10",
                               superscript("-5"))),
                 (value = 1.0e-7, tag = "1e-7",
                  label = rich(rich("rtol", font = RTOL_FONT), " = 10",
                               superscript("-7"))),
                 (value = 1.0e-9, tag = "1e-9",
                  label = rich(rich("rtol", font = RTOL_FONT), " = 10",
                               superscript("-9"))))
const STUDIES = ((name = "richards_celia_haverkamp_mixed", final_time = 360.0),
                 (name = "richards_celia_haverkamp_pressure_head_pressure_head_transfer",
                  final_time = 360.0),
                 (name = "richards_celia_haverkamp_pressure_head_water_content_transfer",
                  final_time = 360.0),
                 (name = "richards_celia_new_mexico_mixed", final_time = 86_400.0),
                 (name = "richards_celia_new_mexico_pressure_head_pressure_head_transfer",
                  final_time = 86_400.0),
                 (name = "richards_celia_new_mexico_pressure_head_water_content_transfer",
                  final_time = 86_400.0))

data_filename(study_name, tolerance) = "$(study_name)_rtol$(tolerance.tag).dat"

function collect_tables(directories)
    files = Dict{String, String}()
    for directory in directories
        data_directory = joinpath(abspath(directory), "data")
        for path in sort(readdir(data_directory; join = true))
            if !endswith(path, ".dat")
                continue
            end
            name = basename(path)
            if haskey(files, name)
                error("Duplicate conservation table: $name")
            end
            files[name] = path
        end
    end

    tables = []
    for study in STUDIES
        paths = [get(files, data_filename(study.name, tolerance), nothing)
                 for tolerance in RELTOLS]
        present = .!isnothing.(paths)
        if any(present) && !all(present)
            error("Incomplete tolerance set for $(study.name)")
        end
        if all(present)
            push!(tables, (; study..., paths = String.(paths)))
        end
    end
    if isempty(tables)
        error("No complete conservation studies found")
    end
    expected = Set(path for table in tables for path in basename.(table.paths))
    unexpected = setdiff(Set(keys(files)), expected)
    if !isempty(unexpected)
        error("Unexpected conservation tables: $(join(sort!(collect(unexpected)), ", "))")
    end

    return tables
end

function shared_bias_limits(tables)
    positive_biases = Float64[]
    for table in tables, path in table.paths
        _, biases = mass_bias_history(path; time_column = "time_s",
                                      mass_balance_column = "mass_bias_m")
        append!(positive_biases, filter(value -> isfinite(value) && value > 0,
                                        abs.(biases)))
    end
    if isempty(positive_biases)
        error("Mass-bias histories contain no positive finite magnitudes")
    end

    lower_exponent = floor(Int, log10(minimum(positive_biases)))
    upper_exponent = ceil(Int, log10(maximum(positive_biases)))
    if lower_exponent == upper_exponent
        upper_exponent += 1
    end
    limits = (10.0^lower_exponent, 10.0^upper_exponent)
    tick_exponents = filter(iseven, lower_exponent:upper_exponent)
    ticks = 10.0 .^ tick_exponents
    tick_labels = [LaTeXString("10^{$exponent}")
                   for exponent in tick_exponents]
    return (; limits, ticks = (ticks, tick_labels))
end

function plot_study(table, output_directory, bias_axis)
    time_ticks = if table.final_time == 360.0
        collect(0.0:60.0:360.0)
    else
        collect(0.0:20_000.0:80_000.0)
    end
    tick_labels = [@sprintf("%.0f", time) for time in time_ticks]
    legend_position = if endswith(table.name, "_mixed")
        (:left, :top)
    else
        (:right, :bottom)
    end
    output_path = joinpath(output_directory, "$(table.name)_mass_bias.pdf")
    fig = plot_mass_bias(table.paths; output_path,
                         time_column = "time_s", mass_balance_column = "mass_bias_m",
                         labels = getproperty.(RELTOLS, :label),
                         xlabel = L"$t$ (s)",
                         ylabel = L"$|\epsilon_{\mathrm{b}}(t)|$ (m)",
                         absolute = true, ynorm = 1.0,
                         legend_position, yscale = log10,
                         xticks = (time_ticks, tick_labels),
                         xlims = (0.0, table.final_time), ylims = bias_axis.limits)
    ax = fig.content[1]
    ax.yticks = bias_axis.ticks
    save(output_path, fig)
    return output_path
end

function plot_conservation(directories...)
    if isempty(directories)
        error("Supply at least one saved run directory")
    end
    tables = collect_tables(directories)
    bias_axis = shared_bias_limits(tables)

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output_directory = joinpath(abspath(first(directories)), "figures", id)
    mkpath(dirname(output_directory))
    mkdir(output_directory)

    paths = [plot_study(table, output_directory, bias_axis) for table in tables]
    println("Figures: $output_directory")
    return paths
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        error("Usage: julia --project=run $(@__FILE__) RUN_DIRECTORY [RUN_DIRECTORY ...]")
    end
    RichardsConservationPlots.plot_conservation(ARGS...)
end
