@doc raw"""
    solution_data_1d(result::ImplicitSolveResult; index = lastindex(result.sol.u), component = 1)
    solution_data_1d(sol::SciMLBase.AbstractODESolution;
                     index = lastindex(sol.u), component = 1, mesh_history = nothing)

Extract a saved one-dimensional solution on Trixi's default DG plotting grid and return
`(; time, x, values, mesh_vertices_x)`. Both element traces are retained at shared
boundaries. The mesh vertices are returned separately, without padding.
The returned arrays do not alias the solution or its recorded mesh history.
This function is available without loading a visualization package.

For `SemidiscretizationImplicit`, components select evolved variables first, followed by
spatial state variables. For the mixed Richards form, component 1 is water content and
component 2 is pressure head; for the pressure-head form, component 1 is pressure head.

An `ImplicitSolveResult` supplies its recorded mesh history. With adaptive meshes, use
`save_mesh_history = true` during the solve to extract earlier saved states on their
original meshes. Without mesh history, extraction uses the current mesh; an earlier
state is valid only if that mesh has not changed.

```julia
snapshot = solution_data_1d(result; index = 1, component = 2)
time, x, values, mesh_vertices_x = snapshot
snapshot = solution_data_1d(result.sol; index = 1, component = 2,
                            mesh_history = result.mesh_history)
```
"""
function solution_data_1d(sol::SciMLBase.AbstractODESolution;
                          index = lastindex(sol.u), component = 1, mesh_history = nothing)
    if !isnothing(mesh_history) && length(mesh_history) != length(sol.u)
        throw(ArgumentError("The mesh history must match the saved states."))
    end
    u_ode = sol.u[index]
    semi = sol.prob.p
    semi_base = semi isa SemidiscretizationImplicit ? semi.semi_base : semi
    nvariables = Trixi.nvariables(semi_base)

    if semi isa SemidiscretizationImplicit
        if component < 1 || component > 2 * nvariables
            throw(ArgumentError("The requested component is outside the solution state."))
        end
        block, local_component = if component <= nvariables
            evolved_variable_block(u_ode, semi), component
        else
            state_variable_block(u_ode, semi), component - nvariables
        end
    else
        block, local_component = u_ode, component
    end

    if isnothing(mesh_history)
        # Extract the curve and mesh vertices from the same PlotData1D snapshot.
        pd = Trixi.PlotData1D(block, semi_base; solution_variables = Trixi.cons2cons)
        x = collect(pd.x)
        values = vec(pd.data[:, local_component])
        mesh_vertices_x = isnothing(pd.mesh_vertices_x) ? Float64[] :
                          collect(pd.mesh_vertices_x)
    else
        # Rebuild the plotting grid from saved values and their original mesh edges.
        mesh_vertices_x = copy(mesh_history[index])
        n_nodes = Trixi.nnodes(semi_base.solver)
        n_elements = length(mesh_vertices_x) - 1
        if length(block) != nvariables * n_nodes * n_elements
            throw(ArgumentError("The saved state does not match its mesh history."))
        end

        nodal_values = reshape(block, nvariables, n_nodes, n_elements)
        unstructured_data = Array{eltype(block)}(undef, n_nodes, n_elements, 1)
        original_nodes = Array{eltype(mesh_vertices_x)}(undef, 1, 2, n_elements)
        for element in 1:n_elements
            original_nodes[1, 1, element] = mesh_vertices_x[element]
            original_nodes[1, 2, element] = mesh_vertices_x[element + 1]
            for node in 1:n_nodes
                unstructured_data[node, element, 1] = nodal_values[local_component, node,
                                                                 element]
            end
        end

        x, data, _ = Trixi.get_data_1d(original_nodes, unstructured_data, nothing, true)
        x = collect(x)
        values = vec(data[:, 1])
    end

    return (; time = sol.t[index], x, values, mesh_vertices_x)
end

function solution_data_1d(result::ImplicitSolveResult; kwargs...)
    return solution_data_1d(result.sol; mesh_history = result.mesh_history, kwargs...)
end
