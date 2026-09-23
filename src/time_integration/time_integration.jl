@doc raw"""
    default_algorithm(ode::SciMLBase.ODEProblem; kwargs...)

Return a recommended OrdinaryDiffEq.jl time integration algorithm for `ode`, suitable for
passing to `SciMLBase.solve`.

The recommended algorithm for `SemidiscretizationImplicit` is `Rodas5P`, an eight-stage,
fifth-order Rosenbrock-Wanner method, and it recomputes the Jacobian after at most one time
step. When `ode` has the sparse Jacobian prototype supplied by [`semidiscretize`](@ref) with
[`SparseJacobian`](@ref), the algorithm uses sparse forward-mode automatic differentiation
with a deterministic analytical colouring and a KLU linear solver. The colouring keeps
the sparse differentiation cache type unchanged when AMR changes a mesh that is large
enough to contain the complete colour palette. Otherwise, the algorithm uses dense
forward-mode automatic differentiation with automatic chunk-size selection and a dense LU
linear solver.

Keyword arguments override these defaults or are forwarded to the `Rodas5P` constructor.
Configure differentiation through `autodiff` and preconditioning through `linsolve`.
Pass `step_limiter` and `stage_limiter` to `solve` or `solve_implicit`.

# References
- Davis, T. A., and Palamadai Natarajan, E. (2010). Algorithm 907: KLU, a direct sparse
  solver for circuit simulation problems. *ACM Transactions on Mathematical Software*,
  37(3), Article 36.
  [DOI: 10.1145/1824801.1824814](https://doi.org/10.1145/1824801.1824814)
- Revels, J., Lubin, M., and Papamarkou, T. (2016). Forward-mode automatic differentiation
  in Julia. *arXiv:1607.07892*.
  [DOI: 10.48550/arXiv.1607.07892](https://doi.org/10.48550/arXiv.1607.07892)
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
function default_algorithm(ode::SciMLBase.ODEProblem{U, T, I, P};
                           autodiff = OrdinaryDiffEqRosenbrock.AutoForwardDiff(),
                           linsolve = ode.f.jac_prototype isa SparseMatrixCSC ?
                                      LinearSolve.KLUFactorization() :
                                      LinearSolve.LUFactorization(),
                           concrete_jac = nothing,
                           max_jac_age = 1,
                           jac_reuse_gamma_tol = 0.03,
                           kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    return OrdinaryDiffEqRosenbrock.Rodas5P(; autodiff, concrete_jac, linsolve, max_jac_age,
                                            jac_reuse_gamma_tol, kwargs...)
end

"""
    pressure_head_out_of_domain(u, semi, t)

Return `true` if any pressure-head degree of freedom in the candidate ODE state `u` is
nonnegative. The call signature matches SciML's `isoutofdomain` predicate and can be
passed directly to [`solve_implicit`](@ref):
```julia
sol = solve_implicit(ode; isoutofdomain = pressure_head_out_of_domain, kwargs...)
```

The domain criterion dispatches on the equations stored by
[`SemidiscretizationImplicit`](@ref). For [`RichardsEquation1D`](@ref), the pressure head
is obtained from the state-variable block, so the criterion applies to both the mixed and
pressure-head forms. Returning `true` rejects the candidate time step; it does not
constrain internal time-integration stages.
"""
@inline function pressure_head_out_of_domain(u, semi::SemidiscretizationImplicit, t)
    _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)
    return pressure_head_out_of_domain(u, semi, t, equations)
end

function variable_block_norm(block, semi::SemidiscretizationImplicit)
    return function (u, t)
        if u isa Number
            return Trixi.ode_norm(u, t)
        end
        variables = block(u, semi)
        return Trixi.ode_norm(variables, t)
    end
end

"""
    default_stepsize_controller(algorithm, ode)

Return HydroTrixi.jl's default adaptive step-size controller for `algorithm` and `ode`.

For `Rodas5P`, return the PI controller used by [`solve_implicit`](@ref). For other
algorithms, return `nothing` so that OrdinaryDiffEq.jl selects the algorithm's default
controller.
"""
default_stepsize_controller(algorithm, ode) = nothing

function default_stepsize_controller(algorithm::OrdinaryDiffEqRosenbrock.Rodas5P,
                                     ode)
    controller_type = typeof(float(first(ode.tspan)))
    return OrdinaryDiffEqCore.PIController(controller_type, algorithm;
                                           # Current- and previous-error exponents
                                           beta1 = 0.14,
                                           beta2 = 0.08,
                                           # Minimum shrink and maximum growth factors
                                           qmin = 0.2,
                                           qmax = 10.0,
                                           # Allow larger growth after the first step
                                           qmax_first_step = 1.0e4,
                                           # Safety factor
                                           gamma = 0.9,
                                           # OrdinaryDiffEq's deadband, which holds the
                                           # time step fixed when the controller gives
                                           # a time-step divisor between 1.0 and 1.2
                                           qsteady_min = 1.0,
                                           qsteady_max = 1.2,
                                           # Previous-error initialization and floor
                                           qoldinit = 1.0e-4)
end

@doc raw"""
    solve_implicit(ode, algorithm=default_algorithm(ode);
                   dt, adaptive=true, abstol=1.0e-11, reltol=1.0e-7,
                   error_control_block=evolved_variable_block,
                   internalnorm=variable_block_norm(error_control_block, ode.p),
                   error_control_mapping=nothing,
                   kwargs...)

Solve `ode` with HydroTrixi.jl's implicit time integration defaults. The initial time step
`dt` is required. The absolute and relative tolerances default to `1.0e-11` and `1.0e-7`,
respectively.

The adaptive defaults use a PI controller configured for the fifth-order
[`default_algorithm`](@ref), with coefficients ``\beta_1=0.14`` and ``\beta_2=0.08``,
safety factor `0.9`, maximum growth factor `10`, maximum shrink factor `0.2`, and initial
and minimum stored previous error `1.0e-4`. The maximum growth factor is `1.0e4` for the
first step-size proposal and `10` thereafter. OrdinaryDiffEq.jl's default steady-step
deadband holds the time step fixed when the controller proposes a time-step divisor
between `1.0` and `1.2`. These parameters are specified explicitly by
[`default_stepsize_controller`](@ref).

The default `error_control_block = evolved_variable_block` restricts adaptive error
control to the stored evolved-variable block and excludes passive diagnostic variables.
Thus, it selects stored water content for [`MixedForm`](@ref) and stored pressure head
for [`PressureHeadForm`](@ref). Use
`error_control_block = state_variable_block` to select the stored state-variable block
instead. The selected block also defines the default internal norm passed to
OrdinaryDiffEq; pass `internalnorm` explicitly to override it.

With an explicit function `error_control_mapping(value, equations)`, apply that function
pointwise to entries selected by `error_control_block` in the candidate, embedded, and
previous solutions before forming and scaling the error. For example, control water
content reconstructed from stored pressure head with:
```julia
sol = solve_implicit(ode; dt = 1.0e-2,
                     error_control_block = state_variable_block,
                     error_control_mapping = water_content,
                     abstol = 1.0e-8, reltol = 1.0e-5)
```
The tolerances then apply to the mapped quantities. Scalar tolerances are broadcast;
array `abstol` must follow the full ODE layout, and entries corresponding to the selected
block are used. `reltol` must be scalar because it is also used by the Rosenbrock linear
solves.

The mapped residual is reduced with `Trixi.ode_norm` (an MPI-aware RMS norm).
The block norm remains available to OrdinaryDiffEq for its original estimator and other
solver operations; it is not applied to the mapped residual. A controller wrapper
replaces the original error estimate before step-size selection and acceptance. The
original error scaling and norm are still computed; stages, Jacobians, and linear solves
are not repeated. `error_control_mapping` has no effect when `adaptive = false`.

Mapped error control supports the in-place implementations of `Rodas4`, `Rodas42`,
`Rodas4P`, `Rodas4P2`, `Rodas5`, `Rodas5P`, `Rodas5Pe`, and `Rodas6P` with the default
step limiter. It wraps an OrdinaryDiffEq controller, such as
`OrdinaryDiffEqCore.PIController(algorithm)`; if `controller` is `nothing`, it creates
that PI controller. Unsupported algorithms are rejected when mapping is enabled.

Without mapping, HydroTrixi.jl's default controller is used only with `Rodas5P`.
For any other integration algorithm, OrdinaryDiffEq.jl selects its default controller
unless `controller` is passed explicitly.
"""
function solve_implicit(ode::SciMLBase.ODEProblem{U, T, I, P},
                        algorithm = default_algorithm(ode);
                        dt,
                        adaptive = true,
                        abstol = 1.0e-11,
                        reltol = 1.0e-7,
                        controller = default_stepsize_controller(algorithm, ode),
                        dtmin = zero(last(ode.tspan) - first(ode.tspan)),
                        dtmax = last(ode.tspan) - first(ode.tspan),
                        force_dtmin = false,
                        failfactor = 2, # not used by Rodas5P (a linearly implicit method)
                        maxiters = typemax(Int),
                        error_control_block = evolved_variable_block,
                        internalnorm = variable_block_norm(error_control_block, ode.p),
                        error_control_mapping = nothing,
                        save_everystep = false,
                        unstable_check = Trixi.mpi_isparallel() ?
                                         Trixi.ode_unstable_check :
                                         DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
                        kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    if adaptive && error_control_mapping !== nothing
        controller = StepsizeControllerMappedError(controller, algorithm, ode,
                                                    error_control_block,
                                                    error_control_mapping; reltol)
    end
    common_options = (; dt, adaptive, dtmin, dtmax, force_dtmin, failfactor,
                      maxiters, internalnorm, save_everystep, unstable_check)
    adaptive_options = adaptive ? (; abstol, reltol, controller) : (;)
    options = merge(common_options, adaptive_options, (; kwargs...))
    return SciMLBase.solve(ode, algorithm; options...)
end
