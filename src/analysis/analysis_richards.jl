# Total water storage is represented differently by each implicit form
@inline function water_content_integral(u_ode::AbstractVector,
                                        semi::SemidiscretizationImplicit,
                                        ::TemporalOperatorConstitutive)
    return first(evolved_variables_integral(u_ode, semi))
end

@inline function water_content_integral(u_ode::AbstractVector,
                                        semi::SemidiscretizationImplicit,
                                        ::TemporalOperatorCapacity)
    return Trixi.integrate(water_content, u_ode, semi; normalize = false)
end

@inline function water_content_integral(u, semi::SemidiscretizationImplicit,
                                        ::TemporalOperatorConstitutive)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return first(Trixi.integrate(Trixi.cons2cons, u, mesh, equations, solver, cache;
                                 normalize = false))
end

@inline function water_content_integral(u, semi::SemidiscretizationImplicit,
                                        ::TemporalOperatorCapacity)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return Trixi.integrate(water_content, u, mesh, equations, solver, cache;
                           normalize = false)
end

@inline function Trixi.analyze(::typeof(water_content), du::AbstractVector,
                               u_ode::AbstractVector, t, semi::SemidiscretizationImplicit)
    return water_content_integral(u_ode, semi, semi.operator_temporal)
end

@inline function Trixi.analyze(::typeof(water_content), du, u, t,
                               semi::SemidiscretizationImplicit)
    return water_content_integral(u, semi, semi.operator_temporal)
end

Trixi.pretty_form_ascii(::typeof(water_content)) = "water_content"
Trixi.pretty_form_utf(::typeof(water_content)) = "∫θ"

@doc raw"""
    water_content_timederivative

Analysis integral for the time derivative of the total water content,
``\mathrm{d} \int \vartheta(\psi) \,\mathrm{d}z / \mathrm{d}t``.

For the mixed formulation, this integrates the derivative of the evolved variable
``\theta``. For the pressure-head formulation, it integrates
``C(\psi) \partial_t\psi`` using [`water_capacity`](@ref). The two expressions agree when
the constitutive constraint ``\theta = \vartheta(\psi)`` is satisfied.
"""
function water_content_timederivative end

# Storage derivatives use evolved water content or capacity-weighted pressure head rates
@inline function water_content_timederivative_integral(du_ode::AbstractVector,
                                                       u_ode::AbstractVector,
                                                       semi::SemidiscretizationImplicit,
                                                       ::TemporalOperatorConstitutive)
    return first(Trixi.integrate(Trixi.cons2cons, du_ode, semi; normalize = false))
end

@inline function water_content_timederivative_integral(du_ode::AbstractVector,
                                                       u_ode::AbstractVector,
                                                       semi::SemidiscretizationImplicit,
                                                       operator_temporal::TemporalOperatorCapacity)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
    du = Trixi.wrap_array(du_ode, mesh, equations, solver, cache)
    return water_content_timederivative_integral(du, u, semi, operator_temporal)
end

@inline function water_content_timederivative_integral(du, u,
                                                       semi::SemidiscretizationImplicit,
                                                       ::TemporalOperatorConstitutive)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return first(Trixi.integrate(Trixi.cons2cons, du, mesh, equations, solver, cache;
                                 normalize = false))
end

@inline function water_content_timederivative_integral(du, u,
                                                       semi::SemidiscretizationImplicit,
                                                       ::TemporalOperatorCapacity)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return Trixi.integrate_via_indices(u, mesh, equations, solver, cache, du;
                                       normalize = false) do u, i, element, equations,
                                                             solver, du
        u_node = Trixi.get_node_vars(u, equations, solver, i, element)
        du_node = Trixi.get_node_vars(du, equations, solver, i, element)
        return water_capacity(u_node, equations) * pressure_head(du_node)
    end
end

@inline function Trixi.analyze(::typeof(water_content_timederivative),
                               du_ode::AbstractVector, u_ode::AbstractVector, t,
                               semi::SemidiscretizationImplicit)
    return water_content_timederivative_integral(du_ode, u_ode, semi,
                                                 semi.operator_temporal)
end

@inline function Trixi.analyze(::typeof(water_content_timederivative), du, u, t,
                               semi::SemidiscretizationImplicit)
    return water_content_timederivative_integral(du, u, semi, semi.operator_temporal)
end

Trixi.pretty_form_ascii(::typeof(water_content_timederivative)) = "water_content_t"
Trixi.pretty_form_utf(::typeof(water_content_timederivative)) = "d/dt ∫θ"

# Mass bias combines passive boundary fluxes with recorded initial storage
const MASS_BIAS_ANALYSIS_CONTEXT = :HydroTrixi_mass_bias_analysis_context
const MASS_BIAS_INITIAL_WATER_CONTENT = IdDict{Any, Any}()

@doc raw"""
    mass_bias
    mass_bias(u_ode, semi::SemidiscretizationImplicit)
    mass_bias(u_ode, semi::SemidiscretizationImplicit, initial_water_content)

Return the mass balance bias
```math
\epsilon_B(t_M) =
\int_{t_0}^{t_M}
\left(\hat{f}_K(t) - \hat{f}_0(t)\right)\,\mathrm{d}t
- \left(M_h(t_M) - M_h(t_0)\right),
```
where ``\hat{f}_0`` and ``\hat{f}_K`` are the numerical fluxes at the soil surface and
bottom of the column, respectively, and ``M_h`` is the quadrature-based discrete water
mass.

This solver flux output method diagnostic requires
[`PassiveVariablesBoundaryFlux1D`](@ref), which advances the integrated numerical
boundary fluxes as passive variables using the same time integrator as the physical
state. The initial total water content is recorded when the implicit coefficients are
initialized and is used by the analysis callback.
"""
function mass_bias end

function record_mass_bias_initial_storage!(u_ode::AbstractVector,
                                           semi::SemidiscretizationImplicit)
    MASS_BIAS_INITIAL_WATER_CONTENT[semi] = water_content_integral(u_ode, semi,
                                                                   semi.operator_temporal)
    return nothing
end

function mass_bias_initial_water_content(semi::SemidiscretizationImplicit)
    return MASS_BIAS_INITIAL_WATER_CONTENT[semi]
end

function mass_bias(u_ode::AbstractVector, semi::SemidiscretizationImplicit,
                   initial_water_content)
    boundary_fluxes = boundary_flux_integrals(u_ode, semi)
    storage_change = water_content_integral(u_ode, semi, semi.operator_temporal) -
                     initial_water_content
    return -boundary_fluxes.x_neg + boundary_fluxes.x_pos - storage_change
end

function mass_bias(u_ode::AbstractVector, semi::SemidiscretizationImplicit)
    return mass_bias(u_ode, semi, mass_bias_initial_water_content(semi))
end

function mass_bias_initial_water_content(sol, semi::SemidiscretizationImplicit,
                                         initial_water_content)
    !isnothing(initial_water_content) && return initial_water_content

    recorded_water_content = get(MASS_BIAS_INITIAL_WATER_CONTENT, semi, nothing)
    !isnothing(recorded_water_content) && return recorded_water_content

    if first(sol.t) == first(sol.prob.tspan)
        return water_content_integral(first(sol.u), semi, semi.operator_temporal)
    end
end

@doc raw"""
    mass_bias_history(sol; initial_water_content = nothing)
    mass_bias_history(analysis_path::AbstractString;
                      time_column = "time", mass_bias_column = "mass_bias")

Return the saved times and corresponding [`mass_bias`](@ref) values for `sol`.

The solution must use a [`SemidiscretizationImplicit`](@ref) with
[`PassiveVariablesBoundaryFlux1D`](@ref). If no initial water content has been recorded
for the semidiscretization and the first saved state is not at the initial time, pass the
initial total water content explicitly with `initial_water_content`.

When `analysis_path` is provided, read the time and mass-bias columns from a Trixi.jl
analysis file written by `AnalysisCallback(save_analysis = true)`.
"""
function mass_bias_history(sol; initial_water_content = nothing)
    semi = sol.prob.p
    initial_storage = mass_bias_initial_water_content(sol, semi, initial_water_content)
    biases = [mass_bias(u_ode, semi, initial_storage) for u_ode in sol.u]

    return collect(sol.t), biases
end

function mass_bias_history(analysis_path::AbstractString; time_column = "time",
                           mass_bias_column = "mass_bias")
    times = Float64[]
    biases = Float64[]

    open(analysis_path, "r") do io
        # Read the Trixi.jl analysis header to locate scalar output columns
        header = nothing
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue
            header = stripped_line
            break
        end

        header_columns = split(strip(header[2:end]))
        column_indices = Dict(column => index for (index, column) in pairs(header_columns))
        time_index = column_indices[time_column]
        mass_bias_index = column_indices[mass_bias_column]

        # Parse the scalar time history from the selected columns
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue
            startswith(stripped_line, "#") && continue

            values = split(stripped_line)
            push!(times, parse(Float64, values[time_index]))
            push!(biases, parse(Float64, values[mass_bias_index]))
        end
    end

    if isempty(times)
        throw(ArgumentError("`analysis_path` does not contain mass-bias samples."))
    end

    return times, biases
end

# AnalysisCallback supplies the full ODE state needed for passive diagnostics
function (analysis_callback::Trixi.AnalysisCallback)(io, du, u, u_ode, t,
                                                     semi::SemidiscretizationImplicit)
    return task_local_storage(MASS_BIAS_ANALYSIS_CONTEXT, (; u_ode)) do
        invoke(analysis_callback, Tuple{Any, Any, Any, Any, Any, Any},
               io, du, u, u_ode, t, semi)
    end
end

@inline function Trixi.analyze(::typeof(mass_bias), du, u, t,
                               semi::SemidiscretizationImplicit)
    state = task_local_storage(MASS_BIAS_ANALYSIS_CONTEXT)
    return mass_bias(state.u_ode, semi)
end

Trixi.pretty_form_ascii(::typeof(mass_bias)) = "mass_bias"
Trixi.pretty_form_utf(::typeof(mass_bias)) = "ε_b"
