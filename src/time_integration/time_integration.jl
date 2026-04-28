@doc raw"""
    default_algorithm(semi)

Return a recommended `OrdinaryDiffEq` time-integration algorithm for the
semidiscretization `semi`, suitable for passing to `SciMLBase.solve`.

The chosen algorithm depends on the type of `semi`:

- For a generic `Trixi.AbstractSemidiscretization`, returns `Tsit5()`, a 5(4)
  adaptive explicit Runge-Kutta method that is a sensible default for
  non-stiff parabolic problems.
- For a [`SemidiscretizationConstitutive`](@ref), returns
  `Rodas5P(autodiff = AutoFiniteDiff())`, a stiffly accurate linearly
  implicit Rosenbrock method.
"""
default_algorithm(::Trixi.AbstractSemidiscretization) = OrdinaryDiffEq.Tsit5()

function default_algorithm(::SemidiscretizationConstitutive)
    OrdinaryDiffEq.Rodas5P(autodiff = OrdinaryDiffEq.OrdinaryDiffEqDifferentiation.AutoFiniteDiff())
end