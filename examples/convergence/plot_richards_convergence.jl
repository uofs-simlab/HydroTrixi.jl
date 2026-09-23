# Replot one or more saved result directories; never runs a simulation.

module RichardsConvergencePlots

using CairoMakie, HydroTrixi, LaTeXStrings
using Dates, Printf

const ERROR_LIMITS = (5.0e-13, 1.0e-3) # pressure-head error in metres
const SPATIAL_X_LIMITS = (1.5e-4, 5.0e-3)
const SPATIAL_X_TICKS = [0.0002, 0.0005, 0.001, 0.002, 0.004]
const SPATIAL_TRIANGLE_GAP = 2.5
const TEMPORAL_TRIANGLE_GAP = 3.0
const FIGURE_SIZE = (500, 350)

function draw_reference_triangle!(ax, groups, order; position = :below, gap_factor = 1.5)
    fine_x, coarse_x = first(groups).x[1:2]
    ratio = coarse_x / fine_x

    # Use both endpoint errors to place the triangle close to the curves.
    bound = position === :below ? minimum : maximum
    combine = position === :below ? min : max
    coarse_error = bound(errors[2] for group in groups for errors in group.errors)
    fine_error = bound(errors[1] for group in groups for errors in group.errors)
    reference_error = combine(coarse_error, fine_error * ratio^order)

    if position === :below
        HydroTrixi.plot_bottom_triangle!(ax, coarse_x, fine_x, reference_error, order;
                                        triangle_slope = :positive, gap_factor,
                                        trianglefontsize = 15)
    else
        # Upside-down triangle: horizontal edge on top, vertical edge on the left.
        upper = reference_error * gap_factor
        lower = upper / ratio^order
        lines!(ax, [fine_x, coarse_x, fine_x, fine_x], [upper, upper, lower, upper];
               color = :black)
        label_x = 10^((2 * log10(fine_x) + log10(coarse_x)) / 3)
        label_y = 10^((2 * log10(upper) + log10(lower)) / 3)
        text!(ax, label_x, label_y; text = "$order:1", align = (:center, :center),
              color = :black, fontsize = 15, font = HydroTrixi.DEFAULT_PLOT_FONT)
    end
    return nothing
end

function read_table(path)
    lines = readlines(path)
    names = Tuple(Symbol.(split(first(lines))))
    return [NamedTuple{names}(Tuple(names[i] in (:form, :retcode, :status) ? s : parse(Float64, s)
                                   for (i, s) in enumerate(split(line))))
            for line in lines[2:end] if !isempty(strip(line))]
end

function plot_study(path, rows, temporal, N, output_dir)
    xfield = temporal ? :dt_s : :delta_z_m
    groups = map((("mixed", 2, :rect, 11.0), ("pressure_head", 1, :circle, 7.0))) do (name, color, marker, markersize)
        # Ascending spacing/time step places the finest level first.
        selected = sort(filter(r -> r.form == name, rows); by = r -> getproperty(r, xfield))
        labels = name == "mixed" ? (L"$L^2$, mixed form", L"$L^\infty$, mixed form") :
                                  (L"$L^2$, pressure-head form", L"$L^\infty$, pressure-head form")
        (; x = [getproperty(r, xfield) for r in selected],
           errors = ([r.l2_m for r in selected], [r.linf_m for r in selected]),
           labels, color, marker, markersize)
    end

    ticks = first(groups).x
    tick_spec = temporal ? (ticks, [@sprintf("%.3g", x) for x in ticks]) :
                          (SPATIAL_X_TICKS, string.(SPATIAL_X_TICKS))
    output_path = joinpath(output_dir, replace(basename(path), ".dat" => ".pdf"))
    fig = plot_convergence_1d(groups; output_path, size = FIGURE_SIZE,
        xlabel = temporal ? L"$\Delta t$ (s)" : L"$\Delta z$ (m)",
        ylabel = "Pressure-head error (m)",
        xticks = tick_spec, xlims = temporal ? nothing : SPATIAL_X_LIMITS,
        ylims = ERROR_LIMITS, legend_position = (:left, :top))
    ax = fig.content[1]

    draw_reference_triangle!(ax, groups, temporal ? 4 : N + 1;
        gap_factor = temporal ? TEMPORAL_TRIANGLE_GAP : SPATIAL_TRIANGLE_GAP)
    if temporal
        draw_reference_triangle!(ax, groups, 5; position = :above,
                                 gap_factor = TEMPORAL_TRIANGLE_GAP)
    end

    # Label each decade on the error axis.
    error_exponents = ceil(Int, log10(first(ERROR_LIMITS))):floor(Int, log10(last(ERROR_LIMITS)))
    ax.yticks = (10.0 .^ error_exponents, [LaTeXString("10^{$p}") for p in error_exponents])

    save(output_path, fig)
    return output_path
end

function plot_convergence(directories...)
    if isempty(directories)
        error("Supply at least one saved run directory")
    end
    tables = []
    for directory in directories
        data = joinpath(abspath(directory), "data")
        for path in sort(readdir(data; join = true))
            if !endswith(path, ".dat")
                continue
            end
            rows = read_table(path)
            temporal = occursin("_time_", basename(path))
            N = Int(first(rows).N)
            push!(tables, (; path, rows, temporal, N))
        end
    end
    if isempty(tables) || !allunique(basename(t.path) for t in tables)
        error("Missing or duplicate studies")
    end

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output_dir = joinpath(abspath(first(directories)), "figures", id)
    mkpath(dirname(output_dir))
    mkdir(output_dir)

    paths = [plot_study(t.path, t.rows, t.temporal, t.N, output_dir) for t in tables]
    println("Figures: $output_dir")
    return paths
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        error("Usage: julia --project=run $(@__FILE__) RUN_DIRECTORY [RUN_DIRECTORY ...]")
    end
    RichardsConvergencePlots.plot_convergence(ARGS...)
end
