@doc raw"""
    default_algorithm(semi)
    default_algorithm(semi, jacobian)

Return a recommended `OrdinaryDiffEq` time-integration algorithm for the
semidiscretization `semi`, suitable for passing to `SciMLBase.solve`. When `jacobian` is
provided, select linear solver settings compatible with that Jacobian strategy.

The chosen algorithm depends on the type of `semi`:

- For a generic `Trixi.AbstractSemidiscretization`, returns `Tsit5()`, a 5(4) adaptive
  explicit Runge-Kutta method that is a sensible default for non-stiff problems.
- For a [`SemidiscretizationImplicit`](@ref) using [`DefaultJacobian`](@ref), returns
  `Rodas5P(autodiff = AutoFiniteDiff())`, an eight-stage, fifth-order Rosenbrock-Wanner
  method with an embedded fourth-order approximation for adaptive time stepping.
- For a [`SemidiscretizationImplicit`](@ref) using [`SparseJacobian`](@ref), additionally
  configures `KLUFactorization()` as the sparse linear solver.
"""
default_algorithm(::Trixi.AbstractSemidiscretization) = OrdinaryDiffEq.Tsit5()

function default_algorithm(::SemidiscretizationImplicit)
    autodiff = OrdinaryDiffEq.OrdinaryDiffEqDifferentiation.AutoFiniteDiff()
    return OrdinaryDiffEq.Rodas5P(autodiff = autodiff)
end

# Keep the linear solver compatible with the Jacobian storage strategy
function default_algorithm(semi::SemidiscretizationImplicit, ::DefaultJacobian)
    default_algorithm(semi)
end

function default_algorithm(::SemidiscretizationImplicit, ::SparseJacobian)
    autodiff = OrdinaryDiffEq.OrdinaryDiffEqDifferentiation.AutoFiniteDiff()
    linsolve = LinearSolve.KLUFactorization()
    return OrdinaryDiffEq.Rodas5P(autodiff = autodiff, linsolve = linsolve)
end
