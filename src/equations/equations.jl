"""
    BoundaryConditionDirichletPenalty(boundary_value_function; penalty=1.0)

Dirichlet boundary condition with a configurable boundary penalty for parabolic
diffusion operators. The `boundary_value_function` is called as
`boundary_value_function(x, t, equations)`.

`penalty` may be either a scalar or a function called as
`penalty(x, t, equations)`.
"""
struct BoundaryConditionDirichletPenalty{B, P}
    boundary_value_function::B
    penalty_function::P
end

function BoundaryConditionDirichletPenalty(boundary_value_function; penalty=1.0)
    penalty_function = penalty isa Function ? penalty : (x, t, equations) -> penalty
    return BoundaryConditionDirichletPenalty(boundary_value_function, penalty_function)
end
