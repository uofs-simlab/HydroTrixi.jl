"""
    examples_dir()

Return the absolute path to the package `examples/` directory.
"""
examples_dir() = joinpath(pkgdir(HydroTrixi), "examples")

"""
    compute_eoc(errors; refinement_factor=2.0)

Compute the estimated orders of convergence between successive entries of `errors`.
The first entry of the returned vector is `NaN` because no coarser reference value is
available.
"""
function compute_eoc(errors::AbstractVector{<:Real}; refinement_factor = 2.0)
    values = similar(errors, Float64)
    fill!(values, NaN)

    denom = log(refinement_factor)

    first_i = firstindex(errors)
    last_i = lastindex(errors)
    for i in (first_i + 1):last_i
        values[i] = log(errors[i - 1] / errors[i]) / denom
    end

    return values
end
