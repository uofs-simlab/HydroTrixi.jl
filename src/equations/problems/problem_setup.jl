@doc raw"""
    HydrologicProblem(; equations, initial_condition, boundary_conditions, domain, tspan,
                        source_terms = nothing, state_to_evolved = nothing,
                        evolved_to_state = nothing)

A container for reusable hydrologic PDE problem data. It stores the governing `equations`,
an `initial_condition`, `boundary_conditions`, the spatial `domain`, the time interval
`tspan`, optional `source_terms`, and optional state-to-evolved and evolved-to-state maps
for mixed-form problems.

The `domain` must be a pair `(x_min, x_max)` of coordinate tuples with matching dimension.
The type parameter `NDIMS` records that dimension for dispatch and introspection.
"""
struct HydrologicProblem{NDIMS, Equations, InitialCondition, BoundaryConditions,
                         Domain, Tspan, SourceTerms, StateToEvolved, EvolvedToState}
    equations::Equations
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    domain::Domain
    tspan::Tspan
    source_terms::SourceTerms
    state_to_evolved::StateToEvolved
    evolved_to_state::EvolvedToState
end

function Base.show(io::IO, hydrologic_problem::HydrologicProblem)
    @nospecialize hydrologic_problem # reduce precompilation time

    print(io, "HydrologicProblem(")
    print(io, hydrologic_problem.equations)
    print(io, ", ", hydrologic_problem.initial_condition)
    print(io, ", ", hydrologic_problem.boundary_conditions)
    print(io, ", ", hydrologic_problem.domain)
    print(io, ", ", hydrologic_problem.tspan)
    print(io, ", ", hydrologic_problem.source_terms)
    print(io, ", ", hydrologic_problem.state_to_evolved)
    print(io, ", ", hydrologic_problem.evolved_to_state)
    print(io, ")")
    return nothing
end

function HydrologicProblem(; equations,
                           initial_condition,
                           boundary_conditions,
                           domain,
                           tspan,
                           source_terms = nothing,
                           state_to_evolved = nothing,
                           evolved_to_state = nothing)
    length(domain) == 2 ||
        throw(ArgumentError("Expected `domain` to be a pair of coordinate tuples."))
    lower, upper = domain
    ndims = length(lower)
    length(upper) == ndims ||
        throw(ArgumentError("Expected lower and upper domain bounds to have " *
                            "matching dimension."))

    return HydrologicProblem{ndims, typeof(equations), typeof(initial_condition),
                             typeof(boundary_conditions), typeof(domain), typeof(tspan),
                             typeof(source_terms), typeof(state_to_evolved),
                             typeof(evolved_to_state)}(equations,
                                                       initial_condition,
                                                       boundary_conditions,
                                                       domain,
                                                       tspan,
                                                       source_terms,
                                                       state_to_evolved,
                                                       evolved_to_state)
end

Base.ndims(::HydrologicProblem{NDIMS}) where {NDIMS} = NDIMS

@inline function boundary_condition_summary_name(boundary_condition)
    return nameof(typeof(boundary_condition))
end

function pretty_boundary_name(boundary_name::Symbol)
    boundary_name === :x_neg && return "negative x"
    boundary_name === :x_pos && return "positive x"
    boundary_name === :y_neg && return "negative y"
    boundary_name === :y_pos && return "positive y"
    boundary_name === :z_neg && return "negative z"
    boundary_name === :z_pos && return "positive z"
    return String(boundary_name)
end

function print_boundary_conditions_summary(io::IO, boundary_conditions::NamedTuple)
    Trixi.summary_line(io, "boundary conditions", length(boundary_conditions))
    for (boundary_name, boundary_condition) in pairs(boundary_conditions)
        Trixi.summary_line(Trixi.increment_indent(io),
                           pretty_boundary_name(boundary_name),
                           boundary_condition_summary_name(boundary_condition))
    end
    return nothing
end

function print_boundary_conditions_summary(io::IO, boundary_conditions)
    return Trixi.summary_line(io, "boundary conditions",
                              boundary_condition_summary_name(boundary_conditions))
end

function Base.show(io::IO, ::MIME"text/plain", hydrologic_problem::HydrologicProblem)
    @nospecialize hydrologic_problem # reduce precompilation time

    if get(io, :compact, false)
        show(io, hydrologic_problem)
    else
        Trixi.summary_header(io, "HydrologicProblem")
        Trixi.summary_line(io, "#spatial dimensions", ndims(hydrologic_problem))
        Trixi.summary_line(io, "equations",
                           hydrologic_problem.equations |> typeof |> nameof)
        Trixi.summary_line(io, "initial condition",
                           hydrologic_problem.initial_condition)
        print_boundary_conditions_summary(io, hydrologic_problem.boundary_conditions)
        Trixi.summary_line(io, "domain", hydrologic_problem.domain)
        Trixi.summary_line(io, "time interval", hydrologic_problem.tspan)
        Trixi.summary_line(io, "source terms", hydrologic_problem.source_terms)
        Trixi.summary_line(io, "state to evolved", hydrologic_problem.state_to_evolved)
        Trixi.summary_line(io, "evolved to state", hydrologic_problem.evolved_to_state)
        Trixi.summary_footer(io)
    end
end

abstract type AbstractImplicitForm end

@doc raw"""
    MixedForm()

Select the mixed form of the Richards equation in
[`SemidiscretizationImplicit`](@ref). The evolved variable is water content, and pressure
head is recovered through the water retention relation.
"""
struct MixedForm <: AbstractImplicitForm end

@doc raw"""
    PressureHeadForm()
    PressureHeadForm(; transfer_variables = pressure_head)

Select the pressure-head form of the Richards equation in
[`SemidiscretizationImplicit`](@ref). The evolved state is pressure head, and the temporal
operator divides the spatial residual by the nodal water capacity ``C(\psi)``. The
optional `transfer_variables` map controls which variable is transferred during adaptive
mesh refinement. The default transfers pressure head directly. Use
`transfer_variables = water_content` to transfer water content and reconstruct pressure
head after each mesh update.
"""
struct PressureHeadForm{TransferVariables} <: AbstractImplicitForm
    transfer_variables::TransferVariables
end

function PressureHeadForm(; transfer_variables = pressure_head)
    return PressureHeadForm(transfer_variables)
end

# Problem forms select the temporal operator used by the implicit semidiscretization
function implicit_temporal_operator(::MixedForm, hydrologic_problem, capacity_function)
    state_to_evolved = hydrologic_problem.state_to_evolved
    evolved_to_state = hydrologic_problem.evolved_to_state
    isnothing(state_to_evolved) &&
        throw(ArgumentError("Hydrologic problem does not define `state_to_evolved`."))
    return TemporalOperatorConstitutive(state_to_evolved;
                                        evolved_to_state = evolved_to_state)
end

function implicit_temporal_operator(form::PressureHeadForm, hydrologic_problem,
                                    capacity_function)
    transfer_variables = form.transfer_variables
    transfer_to_state = transfer_to_state_function(transfer_variables,
                                                   hydrologic_problem)

    return TemporalOperatorCapacity(capacity_function;
                                    transfer_variables = transfer_variables,
                                    transfer_to_state = transfer_to_state)
end

# Pressure-head transfer inverses are inferred from the requested variables
function transfer_to_state_function(transfer_variables, hydrologic_problem)
    if transfer_variables === pressure_head || transfer_variables === Trixi.cons2cons
        return amr_transfer_identity
    end

    if transfer_variables === hydrologic_problem.state_to_evolved
        transfer_to_state = hydrologic_problem.evolved_to_state
        isnothing(transfer_to_state) &&
            throw(ArgumentError("AMR transfer requires `evolved_to_state`."))
        return transfer_to_state
    end

    throw(ArgumentError("Pressure-head AMR transfer requires `transfer_variables` " *
                        "to be `pressure_head`, `Trixi.cons2cons`, or the " *
                        "hydrologic problem `state_to_evolved` map."))
end

@doc raw"""
    SemidiscretizationImplicit(mesh, hydrologic_problem, solver; solver_parabolic,
                               passive_variables = NoPassiveVariables(),
                               form = MixedForm(),
                               source_terms = hydrologic_problem.source_terms,
                               capacity_function = water_capacity, kwargs...)

Create an implicit semidiscretization for `hydrologic_problem` using the parabolic
spatial solver and the selected temporal operator. Use `form = MixedForm()` for the
constitutive mixed form, or `form = PressureHeadForm()` for the capacity-weighted
pressure-head form. The optional `source_terms` keyword defaults to the hydrologic
problem source terms, and `passive_variables` appends diagnostic scalars after the
physical state. Extra keyword arguments are forwarded to
`Trixi.SemidiscretizationParabolic`.
"""
function SemidiscretizationImplicit(mesh, hydrologic_problem::HydrologicProblem,
                                    solver; solver_parabolic,
                                    passive_variables = NoPassiveVariables(),
                                    form = MixedForm(),
                                    source_terms = hydrologic_problem.source_terms,
                                    capacity_function = water_capacity, kwargs...)
    boundary_conditions = hydrologic_problem.boundary_conditions

    semi_base = Trixi.SemidiscretizationParabolic(mesh,
                                                  hydrologic_problem.equations,
                                                  hydrologic_problem.initial_condition,
                                                  solver;
                                                  boundary_conditions = boundary_conditions,
                                                  source_terms = source_terms,
                                                  solver_parabolic = solver_parabolic,
                                                  kwargs...)
    operator_temporal = implicit_temporal_operator(form, hydrologic_problem,
                                                   capacity_function)

    return SemidiscretizationImplicit(semi_base, operator_temporal, passive_variables)
end
