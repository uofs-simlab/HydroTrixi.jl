@doc raw"""
    HydrologicProblem(; equations, initial_condition, boundary_conditions, domain, tspan,
                        state_to_evolved = nothing, evolved_to_state = nothing)

A container for reusable hydrologic PDE problem data. It stores the governing `equations`,
an `initial_condition`, `boundary_conditions`, the spatial `domain`, the time interval
`tspan`, and optional state-to-evolved and evolved-to-state maps for mixed-form problems.

The `domain` must be a pair `(x_min, x_max)` of coordinate tuples with matching dimension.
The type parameter `NDIMS` records that dimension for dispatch and introspection.
"""
struct HydrologicProblem{NDIMS, Equations, InitialCondition, BoundaryConditions,
                         Domain, Tspan, StateToEvolved, EvolvedToState}
    equations::Equations
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    domain::Domain
    tspan::Tspan
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
                             typeof(state_to_evolved),
                             typeof(evolved_to_state)}(equations,
                                                       initial_condition,
                                                       boundary_conditions,
                                                       domain,
                                                       tspan,
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
        Trixi.summary_line(io, "state to evolved", hydrologic_problem.state_to_evolved)
        Trixi.summary_line(io, "evolved to state", hydrologic_problem.evolved_to_state)
        Trixi.summary_footer(io)
    end
end

function SemidiscretizationImplicit(mesh, hydrologic_problem::HydrologicProblem, solver;
                                    solver_parabolic,
                                    kwargs...)
    state_to_evolved = hydrologic_problem.state_to_evolved
    evolved_to_state = hydrologic_problem.evolved_to_state
    boundary_conditions = hydrologic_problem.boundary_conditions
    isnothing(state_to_evolved) &&
        throw(ArgumentError("Hydrologic problem does not define `state_to_evolved`."))

    semi_base = Trixi.SemidiscretizationParabolic(mesh,
                                                  hydrologic_problem.equations,
                                                  hydrologic_problem.initial_condition,
                                                  solver;
                                                  boundary_conditions = boundary_conditions,
                                                  solver_parabolic = solver_parabolic,
                                                  kwargs...)
    operator_temporal = TemporalOperatorConstitutive(state_to_evolved;
                                                     evolved_to_state = evolved_to_state)
    return SemidiscretizationImplicit(semi_base, operator_temporal)
end
