examples_dir() = joinpath(pkgdir(HydroTrixi), "examples")

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
