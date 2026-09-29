@doc raw"""
    animate_solution_1d(sol; output_path = joinpath(pwd(), "solution_1d.mp4"),
                        component = 1, frame_indices = nothing, kwargs...)

Animate a one-dimensional solution history stored in `sol`, write the animation to
`output_path`, and return `output_path`.

Frames are taken from saved solution states. For adaptive meshes, pass `mesh_history`
recorded during the solve so each state uses its original mesh coordinates.
If `frame_indices` is `nothing`, every frame is rendered; otherwise, only the selected
indices are rendered.
Set `component` to choose the plotted variable. Optional keyword arguments control axis
labels and limits, figure size, fonts, line and marker styling, node markers, element
boundary guides, exact-solution overlays, and the output `framerate`.
The horizontal limits default to the column's left and right boundaries.

`exact_solution`, when supplied, is called as `exact_solution(Trixi.SVector(x), t)`.
"""
function HydroTrixi.animate_solution_1d(sol::SciMLBase.AbstractODESolution;
                                        output_path = joinpath(pwd(), "solution_1d.mp4"),
                                        component = 1, framerate = 24,
                                        frame_indices = nothing, mesh_history = nothing,
                                        kwargs...)
    if !isnothing(mesh_history) && length(mesh_history) != length(sol.t)
        throw(ArgumentError("The mesh history must match the saved states."))
    end
    indices = if isnothing(frame_indices)
        collect(eachindex(sol.t))
    else
        collect(frame_indices)
    end
    if isempty(indices)
        throw(ArgumentError("The animation needs at least one frame."))
    end

    first_idx = first(indices)
    x, y, mesh_vertices_x = solution_frame_data_1d(sol, first_idx, mesh_history;
                                                    component = component)
    solution_plot = initialize_solution_plot_1d(x, y, mesh_vertices_x, sol.t[first_idx];
                                                kwargs...)

    mkpath(dirname(abspath(output_path)))

    record(solution_plot.fig, output_path, indices; framerate = framerate,
           px_per_unit = 1) do i
        x_frame, y_frame, mesh_vertices_x_frame =
            solution_frame_data_1d(sol, i, mesh_history; component = component)
        update_solution_plot_1d!(solution_plot, x_frame, y_frame, mesh_vertices_x_frame,
                                 sol.t[i])
    end

    return output_path
end

function HydroTrixi.animate_solution_1d(result::HydroTrixi.ImplicitSolveResult;
                                        kwargs...)
    return HydroTrixi.animate_solution_1d(result.sol;
                                          mesh_history = result.mesh_history, kwargs...)
end
