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
