# Replot one or more saved conservation-study directories; never runs a simulation.

module RichardsConservationPlots

using CairoMakie, HydroTrixi, LaTeXStrings
using Dates, Printf
using DelimitedFiles

const VISUALIZATION = Base.get_extension(HydroTrixi, :HydroTrixiVisualizationExt)

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
const BENCHMARKS = ((name = "haverkamp", final_time = 360.0,
                     solution_ylims = (-0.65, -0.15)),
                    (name = "new_mexico", final_time = 86_400.0,
                     solution_ylims = (-10.5, -0.5)))
const GROUPS = ((name = "mixed", cases = ((name = "mixed", transfer = :water_content),
                                         (name = "mixed_pressure_head_transfer",
                                          transfer = :pressure_head))),
                (name = "pressure_head",
                 cases = ((name = "pressure_head_pressure_head_transfer",
                           transfer = :pressure_head),
                          (name = "pressure_head_water_content_transfer",
                           transfer = :water_content))))

data_stem(benchmark, case) = "richards_celia_$(benchmark)_$(case)"
data_filename(stem, tolerance) = "$(stem)_rtol$(tolerance.tag).dat"
snapshot_filename(stem, tolerance, stage) =
    "$(stem)_rtol$(tolerance.tag)_$(stage)_solution_mesh.dat"

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
    expected = Set{String}()
    for benchmark in BENCHMARKS, group in GROUPS
        transfer_tables = []
        for case in group.cases
            study_name = data_stem(benchmark.name, case.name)
            paths = [get(files, data_filename(study_name, tolerance), nothing)
                     for tolerance in RELTOLS]
            present = .!isnothing.(paths)
            if any(present) && !all(present)
                error("Incomplete tolerance set for $study_name")
            end
            if all(present)
                for path in paths
                    push!(expected, basename(path))
                end
                push!(transfer_tables, (; stem = study_name, transfer = case.transfer,
                                        paths = String.(paths)))
            end
        end
        if !isempty(transfer_tables)
            push!(tables, (; name = data_stem(benchmark.name, group.name),
                           final_time = benchmark.final_time,
                           solution_ylims = benchmark.solution_ylims, transfer_tables))
        end
    end
    if isempty(tables)
        error("No complete conservation studies found")
    end
    unexpected = setdiff(Set(keys(files)), expected)
    if !isempty(unexpected)
        error("Unexpected conservation tables: $(join(sort!(collect(unexpected)), ", "))")
    end

    return tables
end

function accepted_step_history(path)
    columns = split(strip(open(readline, path))[2:end])
    values = readdlm(path, Float64; comments = true)
    steps = values[:, findfirst(==("accepted_step"), columns)]
    times = values[:, findfirst(==("time_s"), columns)]
    dts = values[:, findfirst(==("dt_s"), columns)]
    accepted = steps .> 0
    return times[accepted], dts[accepted]
end

function plot_time_steps(paths, output_path; labels, colors, linestyles,
                         legend_position, xticks, xlims)
    HydroTrixi.set_serif_tex_theme!()
    fig = Figure(size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE, fontsize = 15)
    ax = VISUALIZATION.solution_axis(fig; xlabel = L"$t$ (s)",
                                     ylabel = L"$\Delta t$ (s)", xticks, xlims)
    for (i, path) in pairs(paths)
        times, dts = accepted_step_history(path)
        VISUALIZATION.plot_series!(ax, times, dts; label = labels[i],
                                   color = colors[i], linestyle = linestyles[i])
    end
    VISUALIZATION.add_legend!(ax; position = legend_position,
                              font = HydroTrixi.DEFAULT_PLOT_FONT,
                              labelsize = 14, show_legend = true)
    save(output_path, fig; px_per_unit = 1)
    return output_path
end

function plot_study(table, output_directory, bias_axis)
    time_ticks = if table.final_time == 360.0
        collect(0.0:60.0:360.0)
    else
        collect(0.0:20_000.0:80_000.0)
    end
    tick_labels = [@sprintf("%.0f", time) for time in time_ticks]
    legend_position = endswith(table.name, "_mixed") ? (:right, :top) :
                       (:right, :bottom)
    paths = String[]
    labels = Any[]
    colors = Any[]
    linestyles = Symbol[]
    palette = Makie.wong_colors()
    single_transfer = length(table.transfer_tables) == 1
    for transfer_table in table.transfer_tables
        linestyle = if transfer_table.transfer === :water_content
            :solid
        else
            :dash
        end
        for (i, path) in pairs(transfer_table.paths)
            label = if single_transfer || transfer_table.transfer === :water_content
                RELTOLS[i].label
            else
                nothing
            end
            push!(paths, path)
            push!(labels, label)
            push!(colors, palette[i])
            push!(linestyles, linestyle)
        end
    end

    bias_path = joinpath(output_directory, "$(table.name)_mass_bias.pdf")
    plot_mass_bias_magnitude(paths; output_path = bias_path, time_column = "time_s",
                             mass_balance_column = "mass_bias_m", labels,
                             colors, linestyles, legend_position,
                             xticks = (time_ticks, tick_labels),
                             yticks = bias_axis.ticks,
                             xlims = (0.0, table.final_time),
                             ylims = bias_axis.limits)
    step_path = joinpath(output_directory, "$(table.name)_time_steps.pdf")
    plot_time_steps(paths, step_path; labels, colors, linestyles,
                    legend_position = (:left, :top),
                    xticks = (time_ticks, tick_labels),
                    xlims = (0.0, table.final_time))
    return [bias_path, step_path]
end

function plot_snapshots(table, output_directory)
    output_paths = String[]
    skipped = String[]
    for transfer_table in table.transfer_tables
        for (i, path) in pairs(transfer_table.paths)
            for stage in (:half, :full)
                filename = snapshot_filename(transfer_table.stem, RELTOLS[i], stage)
                pdf_filename = replace(filename, ".dat" => ".pdf")
                snapshot_path = joinpath(dirname(dirname(path)), "snapshots", filename)
                if !isfile(snapshot_path)
                    push!(skipped, pdf_filename)
                    continue
                end
                values = readdlm(snapshot_path, Float64; comments = true)
                t = values[1, 1]
                x = values[:, 2]
                y = values[:, 3]
                mesh_vertices_x = filter(isfinite, values[:, 4])
                output_path = joinpath(output_directory, pdf_filename)
                VISUALIZATION.save_solution_plot_1d(x, y, mesh_vertices_x, t;
                                                    output_path,
                                                    xlabel = "Distance below surface (m)",
                                                    ylabel = "Pressure head (m)",
                                                    ylims = table.solution_ylims,
                                                    show_element_boundaries = true)
                push!(output_paths, output_path)
            end
        end
    end
    if !isempty(skipped)
        println("Skipped $(length(skipped)) snapshot PDFs without saved tables:")
        foreach(path -> println("  $path"), skipped)
    end
    return output_paths
end

function plot_conservation(directories...)
    if isempty(directories)
        error("Supply at least one saved run directory")
    end
    tables = collect_tables(directories)
    paths = [path for table in tables for transfer_table in table.transfer_tables
             for path in transfer_table.paths]
    bias_axis = mass_bias_magnitude_axis(paths; time_column = "time_s",
                                         mass_balance_column = "mass_bias_m")

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output_directory = joinpath(abspath(first(directories)), "figures", id)
    mkpath(dirname(output_directory))
    mkdir(output_directory)

    paths = String[]
    for table in tables
        append!(paths, plot_study(table, output_directory, bias_axis))
        append!(paths, plot_snapshots(table, output_directory))
    end
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
