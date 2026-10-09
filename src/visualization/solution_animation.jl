@doc raw"""
    animate_solution_1d(sol; output_path = joinpath(pwd(), "solution_1d.mp4"),
                        component = 1, frame_indices = nothing, kwargs...)

Animate a one-dimensional solution history stored in `sol`, write the animation to
`output_path`, and return `output_path`.

Frames are taken from saved solution states. For adaptive meshes, pass `mesh_history`
recorded during the solve so each state uses its original mesh coordinates.
If `frame_indices` is `nothing`, every frame is rendered;
otherwise, only the selected indices are rendered.
Set `component` to choose the plotted variable. Optional keyword arguments control axis
labels and limits, figure size, fonts, line and marker styling, node markers, element
boundary guides, exact-solution overlays, and the output `framerate`.
The horizontal limits default to the domain endpoints.

`exact_solution`, when supplied, is called as `exact_solution(Trixi.SVector(x), t)`.
"""
function HydroTrixi.animate_solution_1d(sol::SciMLBase.AbstractODESolution;
                                        output_path = joinpath(pwd(), "solution_1d.mp4"),
                                        component = 1, framerate = 24,
                                        frame_indices = nothing, mesh_history = nothing,
                                        kwargs...)
    indices = if isnothing(frame_indices)
        collect(eachindex(sol.t))
    else
        collect(frame_indices)
    end
    if isempty(indices)
        throw(ArgumentError("The animation needs at least one frame."))
    end

    first_idx = first(indices)
    data = HydroTrixi.solution_data_1d(sol; index = first_idx, component, mesh_history)
    solution_plot = initialize_solution_plot_1d(data.x, data.values, data.mesh_vertices_x,
                                                data.time; kwargs...)

    mkpath(dirname(abspath(output_path)))

    record(solution_plot.fig, output_path, indices; framerate = framerate,
           px_per_unit = 1) do i
        data = HydroTrixi.solution_data_1d(sol; index = i, component, mesh_history)
        update_solution_plot_1d!(solution_plot, data.x, data.values, data.mesh_vertices_x,
                                 data.time)
    end

    return output_path
end

function HydroTrixi.animate_solution_1d(result::HydroTrixi.ImplicitSolveResult;
                                        kwargs...)
    return HydroTrixi.animate_solution_1d(result.sol;
                                          mesh_history = result.mesh_history, kwargs...)
end
