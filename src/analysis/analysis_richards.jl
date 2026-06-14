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
                               u_ode::AbstractVector, t,
                               semi::SemidiscretizationImplicit)
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
``\mathrm{d} \int \theta \,\mathrm{d}x / \mathrm{d}t``.

For the constitutive mixed form, this integrates the evolved water-content derivative.
For the pressure-head form, this integrates
``C(\psi) \partial \psi / \partial t`` using [`water_capacity`](@ref).
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
                               du_ode::AbstractVector,
                               u_ode::AbstractVector, t,
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
\epsilon_B =
\int_{t_0}^{t_M} \left(q(t, 0) - q(t, z_N)\right)\,\mathrm{d}t
- \int_0^L \left(\theta(t_M, z) - \theta(t_0, z)\right)\,\mathrm{d}z.
```

This diagnostic requires [`PassiveVariablesBoundaryFlux1D`](@ref), which stores the
time-integrated boundary fluxes in the DG surface orientation. The initial total water
content is recorded when the implicit coefficients are initialized and is used by the
analysis callback.
"""
function mass_bias end

function record_mass_bias_initial_storage!(u_ode::AbstractVector,
                                           semi::SemidiscretizationImplicit)
    semi.passive_variables isa PassiveVariablesBoundaryFlux1D || return nothing
    MASS_BIAS_INITIAL_WATER_CONTENT[semi] = water_content_integral(u_ode, semi,
                                                                   semi.operator_temporal)
    return nothing
end

function mass_bias_initial_water_content(semi::SemidiscretizationImplicit)
    initial_water_content = get(MASS_BIAS_INITIAL_WATER_CONTENT, semi, nothing)
    isnothing(initial_water_content) &&
        throw(ArgumentError("`mass_bias` requires initialized " *
                            "`PassiveVariablesBoundaryFlux1D` diagnostics."))
    return initial_water_content
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

    throw(ArgumentError("`mass_bias_history` requires either an initialized " *
                        "`PassiveVariablesBoundaryFlux1D` diagnostic, a saved " *
                        "initial state, or `initial_water_content`."))
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

When `analysis_path` is provided, read the time and mass-bias columns from a Trixi
analysis file written by `AnalysisCallback(save_analysis = true)`.
"""
function mass_bias_history(sol; initial_water_content = nothing)
    semi = sol.prob.p
    semi isa SemidiscretizationImplicit ||
        throw(ArgumentError("`mass_bias_history` requires `SemidiscretizationImplicit`."))
    semi.passive_variables isa PassiveVariablesBoundaryFlux1D ||
        throw(ArgumentError("`mass_bias_history` requires " *
                            "`PassiveVariablesBoundaryFlux1D`."))
    isempty(sol.u) && throw(ArgumentError("`sol` does not contain saved states."))

    initial_storage = mass_bias_initial_water_content(sol, semi, initial_water_content)
    biases = [mass_bias(u_ode, semi, initial_storage) for u_ode in sol.u]

    return collect(sol.t), biases
end

function mass_bias_history(analysis_path::AbstractString;
                           time_column = "time",
                           mass_bias_column = "mass_bias")
    isfile(analysis_path) ||
        throw(ArgumentError("`analysis_path` must refer to an existing file."))

    times = Float64[]
    biases = Float64[]

    open(analysis_path, "r") do io
        # Read the Trixi analysis header to locate scalar output columns
        header = nothing
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue
            startswith(stripped_line, "#") ||
                throw(ArgumentError("`analysis_path` does not contain a Trixi " *
                                    "analysis header."))
            header = stripped_line
            break
        end

        isnothing(header) &&
            throw(ArgumentError("`analysis_path` does not contain a header."))

        header_columns = split(strip(header[2:end]))
        time_index = findfirst(==(time_column), header_columns)
        mass_bias_index = findfirst(==(mass_bias_column), header_columns)
        isnothing(time_index) &&
            throw(ArgumentError("`analysis_path` does not contain column " *
                                "`$(time_column)`."))
        isnothing(mass_bias_index) &&
            throw(ArgumentError("`analysis_path` does not contain column " *
                                "`$(mass_bias_column)`."))

        # Parse the scalar time history from the selected columns
        required_columns = max(time_index, mass_bias_index)
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue
            startswith(stripped_line, "#") && continue

            values = split(stripped_line)
            length(values) >= required_columns ||
                throw(ArgumentError("`analysis_path` has a row with too few columns."))
            push!(times, parse(Float64, values[time_index]))
            push!(biases, parse(Float64, values[mass_bias_index]))
        end
    end

    isempty(times) &&
        throw(ArgumentError("`analysis_path` does not contain mass-bias samples."))

    return times, biases
end

function mass_bias_analysis_state()
    state = get(task_local_storage(), MASS_BIAS_ANALYSIS_CONTEXT, nothing)
    isnothing(state) &&
        throw(ArgumentError("`mass_bias` must be evaluated from `AnalysisCallback`."))
    return state
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
    state = mass_bias_analysis_state()
    return mass_bias(state.u_ode, semi)
end

Trixi.pretty_form_ascii(::typeof(mass_bias)) = "mass_bias"
Trixi.pretty_form_utf(::typeof(mass_bias)) = "ε_b"
