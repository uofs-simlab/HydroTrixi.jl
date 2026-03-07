struct LinearDiffusionEquation1D{RealT <: Real} <: Trixi.AbstractLaplaceDiffusion{1, 1}
    diffusivity::RealT
end

Trixi.varnames(::typeof(Trixi.cons2cons), ::LinearDiffusionEquation1D) = ("scalar",)

@inline function Trixi.flux(u, gradients, orientation::Integer,
                            equations::LinearDiffusionEquation1D)
    dudx, = gradients
    return equations.diffusivity * dudx
end
