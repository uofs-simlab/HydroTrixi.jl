mutable struct SemidiscretizationParabolic{Mesh, Equations, InitialCondition,
                                           BoundaryConditions, SourceTerms,
                                           Solver, ParabolicScheme, Cache} <:
               Trixi.AbstractSemidiscretization
    mesh::Mesh
    equations::Equations
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    source_terms::SourceTerms
    solver::Solver
    parabolic_scheme::ParabolicScheme
    cache::Cache
    performance_counter::Trixi.PerformanceCounter
end

function SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                     parabolic_scheme=Trixi.ViscousFormulationLocalDG(),
                                     source_terms=nothing,
                                     boundary_conditions=Trixi.boundary_condition_periodic,
                                     RealT=real(solver),
                                     uEltype=RealT)
    @assert ndims(mesh) == ndims(equations)

    base_cache = Trixi.create_cache(mesh, equations, solver, RealT, uEltype)
    _boundary_conditions = Trixi.digest_boundary_conditions(boundary_conditions, mesh,
                                                            solver, base_cache)
    Trixi.check_periodicity_mesh_boundary_conditions(mesh, _boundary_conditions)

    viscous_cache = Trixi.create_cache_parabolic(mesh, equations, solver,
                                                 Trixi.nelements(solver, base_cache),
                                                 uEltype)
    cache = (; base=base_cache, viscous=viscous_cache)

    performance_counter = Trixi.PerformanceCounter()

    return SemidiscretizationParabolic{typeof(mesh), typeof(equations),
                                       typeof(initial_condition),
                                       typeof(_boundary_conditions),
                                       typeof(source_terms),
                                       typeof(solver), typeof(parabolic_scheme),
                                       typeof(cache)}(mesh, equations,
                                                      initial_condition,
                                                      _boundary_conditions,
                                                      source_terms,
                                                      solver,
                                                      parabolic_scheme,
                                                      cache,
                                                      performance_counter)
end

@inline Base.ndims(semi::SemidiscretizationParabolic) = ndims(semi.mesh)

@inline Trixi.nvariables(semi::SemidiscretizationParabolic) = Trixi.nvariables(semi.equations)

@inline Base.real(semi::SemidiscretizationParabolic) = real(semi.solver)

@inline function Trixi.mesh_equations_solver_cache(semi::SemidiscretizationParabolic)
    return semi.mesh, semi.equations, semi.solver, semi.cache.base
end

function Trixi.calc_error_norms(func, u_ode, t, analyzer,
                                semi::SemidiscretizationParabolic, cache_analysis)
    mesh = semi.mesh
    equations = semi.equations
    initial_condition = semi.initial_condition
    solver = semi.solver
    cache = semi.cache.base
    u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)

    return Trixi.calc_error_norms(func, u, t, analyzer, mesh, equations,
                                  initial_condition, solver, cache, cache_analysis)
end

function Trixi.compute_coefficients(t, semi::SemidiscretizationParabolic)
    return Trixi.compute_coefficients(semi.initial_condition, t, semi)
end

function Trixi.compute_coefficients!(u_ode, t, semi::SemidiscretizationParabolic)
    return Trixi.compute_coefficients!(u_ode, semi.initial_condition, t, semi)
end

function Trixi.rhs!(du_ode, u_ode, semi::SemidiscretizationParabolic, t)
    mesh = semi.mesh
    equations = semi.equations
    boundary_conditions = semi.boundary_conditions
    source_terms = semi.source_terms
    solver = semi.solver
    parabolic_scheme = semi.parabolic_scheme
    base = semi.cache.base
    viscous = semi.cache.viscous

    u = Trixi.wrap_array(u_ode, mesh, equations, solver, base)
    du = Trixi.wrap_array(du_ode, mesh, equations, solver, base)

    time_start = time_ns()
    Trixi.rhs_parabolic!(du, u, t, mesh, equations, boundary_conditions, source_terms,
                         solver, parabolic_scheme, base, viscous)
    runtime = time_ns() - time_start
    put!(semi.performance_counter, runtime)

    return nothing
end
