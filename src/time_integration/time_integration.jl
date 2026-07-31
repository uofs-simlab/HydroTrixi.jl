@doc raw"""
    default_algorithm(semi)

Return a recommended `OrdinaryDiffEq` time integration algorithm for the
semidiscretization `semi`, suitable for passing to `SciMLBase.solve`.

The chosen algorithm depends on the type of `semi`:

- For a generic `Trixi.AbstractSemidiscretization`, returns `Tsit5()`, a 5(4) adaptive
  explicit Runge-Kutta method that is a sensible default for non-stiff problems.
- For a [`SemidiscretizationImplicit`](@ref), returns
  `Rodas5P(autodiff = AutoForwardDiff())`, an eight-stage, fifth-order Rosenbrock-Wanner
  method with graph-coloured forward-mode automatic differentiation when
  [`semidiscretize`](@ref) supplies the default [`SparseJacobian`](@ref) prototype.

The Jacobian strategy passed to [`semidiscretize`](@ref) configures Jacobian storage. The
time integration algorithm configures the differentiation backend and linear solver.

# References
- Tsitouras, C. (2011). Runge-Kutta pairs of order 5(4) satisfying only the first
  column simplifying assumption. *Computers & Mathematics with Applications*, 62(2),
  770-775.
  [DOI: 10.1016/j.camwa.2011.06.002](https://doi.org/10.1016/j.camwa.2011.06.002)
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
default_algorithm(::Trixi.AbstractSemidiscretization) = OrdinaryDiffEq.Tsit5()

# Select the default implicit Rosenbrock algorithm
function default_algorithm(::SemidiscretizationImplicit)
    autodiff = OrdinaryDiffEq.OrdinaryDiffEqDifferentiation.AutoForwardDiff()
    return OrdinaryDiffEq.Rodas5P(autodiff = autodiff)
end
