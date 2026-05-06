@doc raw"""
    HydrologicProblem(; equations, initial_condition, boundary_conditions, domain, tspan,
                        constitutive_relation = nothing)

A container for reusable hydrologic PDE problem data. It stores the governing `equations`,
an `initial_condition`, `boundary_conditions`, the spatial `domain`, the time interval
`tspan`, and an optional `constitutive_relation` for mixed-form Richards' equation problems.

The `domain` must be a pair `(x_min, x_max)` of coordinate tuples with matching dimension.
The type parameter `NDIMS` records that dimension for dispatch and introspection.
"""
struct HydrologicProblem{NDIMS, Equations, InitialCondition, BoundaryConditions,
                         Domain, Tspan, ConstitutiveRelation}
    equations::Equations
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    domain::Domain
    tspan::Tspan
    constitutive_relation::ConstitutiveRelation
end

function Base.show(io::IO, hydrologic_problem::HydrologicProblem)
    @nospecialize hydrologic_problem # reduce precompilation time

    print(io, "HydrologicProblem(")
    print(io, hydrologic_problem.equations)
    print(io, ", ", hydrologic_problem.initial_condition)
    print(io, ", ", hydrologic_problem.boundary_conditions)
    print(io, ", ", hydrologic_problem.domain)
    print(io, ", ", hydrologic_problem.tspan)
    print(io, ", ", hydrologic_problem.constitutive_relation)
    print(io, ")")
    return nothing
end

function HydrologicProblem(; equations,
                           initial_condition,
                           boundary_conditions,
                           domain,
                           tspan,
                           constitutive_relation = nothing)
    length(domain) == 2 ||
        throw(ArgumentError("Expected `domain` to be a pair of coordinate tuples."))
    lower, upper = domain
    ndims = length(lower)
    length(upper) == ndims ||
        throw(ArgumentError("Expected lower and upper domain bounds to have matching dimension."))

    return HydrologicProblem{ndims, typeof(equations), typeof(initial_condition),
                             typeof(boundary_conditions), typeof(domain), typeof(tspan),
                             typeof(constitutive_relation)}(equations,
                                                            initial_condition,
                                                            boundary_conditions,
                                                            domain,
                                                            tspan,
                                                            constitutive_relation)
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
        Trixi.summary_line(io, "constitutive relation",
                           hydrologic_problem.constitutive_relation)
        Trixi.summary_footer(io)
    end
end

function SemidiscretizationImplicit(mesh, hydrologic_problem::HydrologicProblem, solver;
                                    solver_parabolic,
                                    kwargs...)
    constitutive_relation = hydrologic_problem.constitutive_relation
    isnothing(constitutive_relation) &&
        throw(ArgumentError("Hydrologic problem does not define a constitutive relation."))

    semi_base = Trixi.SemidiscretizationParabolic(mesh,
                                                  hydrologic_problem.equations,
                                                  hydrologic_problem.initial_condition,
                                                  solver;
                                                  boundary_conditions = hydrologic_problem.boundary_conditions,
                                                  solver_parabolic = solver_parabolic,
                                                  kwargs...)
    operator_temporal = TemporalOperatorConstitutive(constitutive_relation)
    return SemidiscretizationImplicit(semi_base, operator_temporal)
end
