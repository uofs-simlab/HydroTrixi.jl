@doc raw"""
    default_algorithm(semi::SemidiscretizationImplicit)

Return a recommended OrdinaryDiffEq.jl time integration algorithm for the implicit
semidiscretization `semi`, suitable for passing to `SciMLBase.solve`.

The recommended algorithm is `Rodas5P()`, an eight-stage, fifth-order Rosenbrock-Wanner
method with graph-coloured forward-mode automatic differentiation when
[`semidiscretize`](@ref) supplies the default [`SparseJacobian`](@ref) prototype. The Jacobian strategy passed to [`semidiscretize`](@ref) configures Jacobian storage. The
time integration algorithm configures the differentiation backend and linear solver.

# References
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
# Select the default implicit Rosenbrock algorithm
function default_algorithm(::SemidiscretizationImplicit)
    return OrdinaryDiffEqRosenbrock.Rodas5P()
end
