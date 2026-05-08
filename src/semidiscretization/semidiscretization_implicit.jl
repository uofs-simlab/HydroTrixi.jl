abstract type AbstractTemporalOperator end

@doc raw"""
    SemidiscretizationImplicit{Semidiscretization, TemporalOperator}

A semidiscretization wrapper that augments a spatial semidiscretization with a temporal
operator, so that the semi-discrete problem is not restricted to the explicit form
``\mathrm{d}\mathbf{u}/\mathrm{d}t = R(\mathbf{u}, t)``. Currently, the only supported
temporal operator is [`TemporalOperatorConstitutive`](@ref), which allows for a map
between the state variable(s) and the evolved variable(s) to be enforced in the system of
equations.
"""
struct SemidiscretizationImplicit{Semidiscretization <: Trixi.AbstractSemidiscretization,
                                  TemporalOperator <: AbstractTemporalOperator} <:
       Trixi.AbstractSemidiscretization
    semi_base::Semidiscretization
    operator_temporal::TemporalOperator
end

function Base.show(io::IO, semi::SemidiscretizationImplicit)
    @nospecialize semi # reduce precompilation time

    print(io, "SemidiscretizationImplicit(")
    print(io, semi.semi_base)
    print(io, ", ", semi.operator_temporal |> typeof |> nameof)
    print(io, ")")
    return nothing
end

@doc raw"""
    TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing)

Temporal operator for a [`SemidiscretizationImplicit`](@ref) that takes the form
```math
\frac{\mathrm{d}\mathbf{u}_\mathrm{evolved}}{\mathrm{d} t} =
R(\mathbf{u}_\mathrm{state}, t), \qquad
0 = \mathbf{u}_\mathrm{evolved} - S(\mathbf{u}_\mathrm{state}),
```
where ``R`` is the operator associated with the spatial discretization `semi_base`, and
``S`` is `state_to_evolved`, which maps the state variable(s) to the evolved variable(s).
The optional `evolved_to_state` inverse is required by AMR to reconstruct the algebraic
state after transferring the evolved block. In the case of the Richards equation, for
example, the state variable is the pressure head, and the evolved variable is the water
content, which are related by the water retention curve (i.e., [`water_content`](@ref)).
HydroTrixi.jl then uses SciML's DAE solvers to integrate the resulting system of
equations based on the following mass-matrix form:
```math
\begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix} \frac{\mathrm{d}}{\mathrm{d} t}
\begin{bmatrix} \mathbf{u}_\mathrm{evolved} \\ \mathbf{u}_\mathrm{state} \end{bmatrix} =
\begin{bmatrix}
R(\mathbf{u}_\mathrm{state}, t) \\
\mathbf{u}_\mathrm{evolved} - S(\mathbf{u}_\mathrm{state})
\end{bmatrix}.
"""
struct TemporalOperatorConstitutive{StateToEvolved, EvolvedToState} <:
       AbstractTemporalOperator
    state_to_evolved::StateToEvolved
    evolved_to_state::EvolvedToState
end

function TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing)
    return TemporalOperatorConstitutive(state_to_evolved, evolved_to_state)
end

print_temporal_operator_summary(io::IO, ::AbstractTemporalOperator) = nothing

function print_temporal_operator_summary(io::IO,
                                         operator_temporal::TemporalOperatorConstitutive)
    Trixi.summary_line(io, "state to evolved", operator_temporal.state_to_evolved)
    return Trixi.summary_line(io, "evolved to state",
                              operator_temporal.evolved_to_state)
end

# Wrapper to drive dispatch based on the temporal operator type on methods that take
# mesh, equations, solver, and cache as separate arguments.
struct CacheImplicit{Cache, TemporalOperator <: AbstractTemporalOperator}
    cache_base::Cache
    operator_temporal::TemporalOperator
end

# Most properties are inherited from the base semidiscretization.
@inline function Base.getproperty(semi::SemidiscretizationImplicit, field::Symbol)
    if field === :performance_counter
        return getproperty(getfield(semi, :semi_base), field)
    end
    return getfield(semi, field)
end

@inline function Base.getproperty(cache::CacheImplicit, field::Symbol)
    if field === :cache_base || field === :operator_temporal
        return getfield(cache, field)
    end
    return getproperty(getfield(cache, :cache_base), field)
end

@inline Base.ndims(semi::SemidiscretizationImplicit) = ndims(semi.semi_base)
@inline Base.real(semi::SemidiscretizationImplicit) = real(semi.semi_base)

function Base.show(io::IO, ::MIME"text/plain", semi::SemidiscretizationImplicit)
    @nospecialize semi # reduce precompilation time

    if get(io, :compact, false)
        show(io, semi)
    else
        semi_base = semi.semi_base

        Trixi.summary_header(io, "SemidiscretizationImplicit")
        Trixi.summary_line(io, "#spatial dimensions", ndims(semi))
        Trixi.summary_line(io, "mesh", semi_base.mesh)
        Trixi.summary_line(io, "equations", semi_base.equations |> typeof |> nameof)
        Trixi.summary_line(io, "initial condition", semi_base.initial_condition)
        print_boundary_conditions_summary(io, semi_base.boundary_conditions)
        Trixi.summary_line(io, "source terms", semi_base.source_terms)
        Trixi.summary_line(io, "solver", semi_base.solver |> typeof |> nameof)
        Trixi.summary_line(io, "parabolic solver",
                           semi_base.solver_parabolic |> typeof |> nameof)
        Trixi.summary_line(io, "temporal operator",
                           semi.operator_temporal |> typeof |> nameof)
        print_temporal_operator_summary(io, semi.operator_temporal)
        Trixi.summary_line(io, "total #DOFs per field", Trixi.ndofsglobal(semi))
        Trixi.summary_footer(io)
    end
end

@inline evolved_variable_view(u_ode) = @view(u_ode[1:(length(u_ode) ÷ 2)])
@inline state_variable_view(u_ode) = @view(u_ode[(length(u_ode) ÷ 2 + 1):end])

# Used by analysis callbacks to ensure they extract the evolved variable in the case of a
# `TemporalOperatorConstitutive`.
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any,
                                                       <:TemporalOperatorConstitutive})
    return invoke(Trixi.wrap_array,
                  Tuple{AbstractVector, Trixi.AbstractMesh, Any, Trixi.DGSEM, Any},
                  evolved_variable_view(u_ode), mesh, equations, dg, cache.cache_base)
end

# Default operator hooks correspond to the standard semidiscrete form `∂_t u = R(u, t)`.
@inline function nvariables_total(::AbstractTemporalOperator,
                                  semi_base::Trixi.AbstractSemidiscretization)
    return Trixi.nvariables(semi_base)
end

@inline check_ode_state(u_ode, ::AbstractTemporalOperator) = nothing

@inline function rhs_implicit!(du_ode, u_ode,
                               ::AbstractTemporalOperator,
                               semi_base::Trixi.AbstractSemidiscretization, t)
    return Trixi.default_rhs(semi_base)(du_ode, u_ode, semi_base, t)
end

@inline function mass_matrix(u_ode, ::AbstractTemporalOperator,
                             semi_base::Trixi.AbstractSemidiscretization)
    return Diagonal(ones(eltype(u_ode), length(u_ode)))
end

@inline function nvariables_total(::TemporalOperatorConstitutive,
                                  semi_base::Trixi.AbstractSemidiscretization)
    return 2 * Trixi.nvariables(semi_base)
end

@inline function check_ode_state(u_ode, ::TemporalOperatorConstitutive)
    iseven(length(u_ode)) ||
        throw(ArgumentError("Expected an even number of implicit degrees of freedom."))
    return nothing
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorConstitutive,
                       semi_base::Trixi.AbstractSemidiscretization, t)
    evolved_variable = evolved_variable_view(u_ode)
    state_variable = state_variable_view(u_ode)
    evolved_variable_rhs = evolved_variable_view(du_ode)
    state_variable_rhs = state_variable_view(du_ode)

    # First block of du_ode: R(u_state, t)
    Trixi.default_rhs(semi_base)(evolved_variable_rhs, state_variable, semi_base, t)
    (; equations) = semi_base

    # Second block of du_ode: u_evolved - state_to_evolved(u_state)
    state_to_evolved = operator_temporal.state_to_evolved
    @inbounds for i in eachindex(state_variable_rhs, evolved_variable, state_variable)
        state_variable_rhs[i] = evolved_variable[i] -
                                state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

function mass_matrix(u_ode, ::TemporalOperatorConstitutive,
                     semi_base::Trixi.AbstractSemidiscretization)
    half = length(u_ode) ÷ 2
    diagonal_entries = zeros(eltype(u_ode), length(u_ode))
    @inbounds diagonal_entries[1:half] .= one(eltype(u_ode))
    return Diagonal(diagonal_entries)
end

@inline function Trixi.mesh_equations_solver_cache(semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return mesh, equations, solver, CacheImplicit(cache, semi.operator_temporal)
end

@inline Trixi.default_rhs(::SemidiscretizationImplicit) = rhs_implicit!

# Number of variables depends on the temporal operator type, so we dispatch on that
@inline function Trixi.nvariables(semi::SemidiscretizationImplicit)
    return nvariables_total(semi.operator_temporal, semi.semi_base)
end

function Trixi.compute_coefficients(t, semi::SemidiscretizationImplicit)
    semi_base = semi.semi_base
    operator_temporal = semi.operator_temporal
    coefficients_state = Trixi.compute_coefficients(t, semi_base)
    coefficients_evolved = similar(coefficients_state)
    equations = semi_base.equations
    state_to_evolved = operator_temporal.state_to_evolved

    @inbounds for i in eachindex(coefficients_evolved, coefficients_state)
        coefficients_evolved[i] = state_to_evolved(coefficients_state[i], equations)
    end

    return vcat(coefficients_evolved, coefficients_state)
end

function Trixi.compute_coefficients!(u_ode, t, semi::SemidiscretizationImplicit)
    semi_base = semi.semi_base
    operator_temporal = semi.operator_temporal
    check_ode_state(u_ode, operator_temporal)

    evolved_variable = evolved_variable_view(u_ode)
    state_variable = state_variable_view(u_ode)
    Trixi.compute_coefficients!(state_variable, t, semi_base)
    equations = semi_base.equations
    state_to_evolved = operator_temporal.state_to_evolved

    @inbounds for i in eachindex(evolved_variable, state_variable)
        evolved_variable[i] = state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

function rhs_implicit!(du_ode, u_ode, semi::SemidiscretizationImplicit, t)
    return rhs_implicit!(du_ode, u_ode, semi.operator_temporal, semi.semi_base, t)
end

function Trixi.semidiscretize(semi::SemidiscretizationImplicit, tspan;
                              reset_threads = true)
    if reset_threads
        Trixi.Polyester.reset_threads!()
    end

    u0_ode = Trixi.compute_coefficients(first(tspan), semi)
    check_ode_state(u0_ode, semi.operator_temporal)

    mass_matrix_implicit = mass_matrix(u0_ode, semi.operator_temporal, semi.semi_base)
    ode_function = SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}(rhs_implicit!;
                                                                         mass_matrix = mass_matrix_implicit)
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(ode_function,
                                                                u0_ode,
                                                                tspan,
                                                                semi)
end
