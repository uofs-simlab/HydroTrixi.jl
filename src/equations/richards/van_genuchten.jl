@muladd begin
#! format: noindent

@doc raw"""
    VanGenuchten(; saturated_hydraulic_conductivity, alpha, n, theta_s, theta_r,
                   m = 1 - 1 / n, pore_connectivity = 1 / 2)

A van Genuchten-Mualem soil-hydraulic model,
written in terms of the effective saturation following Ireson et al. (2023).
For ``\psi < 0``, the water content and hydraulic conductivity are
```math
\vartheta(\psi) = \theta_{\mathrm{r}}+
(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})
\left(1+(\alpha_{\mathrm{V}}|\psi|)^{n_{\mathrm{V}}}\right)^{-m_{\mathrm{V}}},
\qquad
\kappa(\psi) = \kappa_{\mathrm{s}} S_{\mathrm{e}}(\psi)^{\ell_{\mathrm{V}}}
\left(1-\left(1-S_{\mathrm{e}}(\psi)^{1/m_{\mathrm{V}}}\right)^{m_{\mathrm{V}}}\right)^2,
```
where
```math
S_{\mathrm{e}}(\psi) \coloneqq
\frac{\vartheta(\psi)-\theta_{\mathrm{r}}}
{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}.
```
The inverse pressure-head scale ``\alpha_{\mathrm{V}}``, water-retention exponent
``n_{\mathrm{V}}``, derived shape parameter ``m_{\mathrm{V}}``, and pore-connectivity
parameter ``\ell_{\mathrm{V}}`` correspond to `alpha`, `n`, `m`, and `pore_connectivity`,
respectively.
The constructor defaults to ``m_{\mathrm{V}}=1-1/n_{\mathrm{V}}``
but permits an explicit `m`.
The parameters ``\kappa_{\mathrm{s}}``, ``\theta_{\mathrm{s}}``,
and ``\theta_{\mathrm{r}}`` correspond to `saturated_hydraulic_conductivity`,
`theta_s`, and `theta_r`, respectively.
For ``\psi\geq 0``, the model uses the saturated values ``S_{\mathrm{e}}(\psi)=1``,
``\vartheta(\psi)=\theta_{\mathrm{s}}``, and ``\kappa(\psi)=\kappa_{\mathrm{s}}``.

# References
- van Genuchten, M. Th. (1980).
  A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
  *Soil Science Society of America Journal*, 44(5), 892-898.
  [DOI: 10.2136/sssaj1980.03615995004400050002x](https://doi.org/10.2136/sssaj1980.03615995004400050002x)
- Mualem, Y. (1976).
  A new model for predicting the hydraulic conductivity of unsaturated porous media.
  *Water Resources Research*, 12(3), 513-522.
  [DOI: 10.1029/WR012i003p00513](https://doi.org/10.1029/WR012i003p00513)
"""
struct VanGenuchten{RealT}
    saturated_hydraulic_conductivity::RealT
    alpha::RealT
    n::RealT
    m::RealT
    pore_connectivity::RealT
    theta_s::RealT
    theta_r::RealT
end

function VanGenuchten(; saturated_hydraulic_conductivity, alpha, n, theta_s, theta_r,
                      m = 1 - 1 / n, pore_connectivity = 1 / 2)
    RealT = promote_type(typeof(saturated_hydraulic_conductivity), typeof(alpha),
                         typeof(n), typeof(m), typeof(pore_connectivity),
                         typeof(theta_s), typeof(theta_r))
    return VanGenuchten{RealT}(saturated_hydraulic_conductivity, alpha, n, m,
                               pore_connectivity, theta_s, theta_r)
end

@inline function effective_saturation(psi, model::VanGenuchten)
    if psi >= zero(psi)
        return one(psi)
    end

    return (one(psi) + (model.alpha * abs(psi))^model.n)^(-model.m)
end

@inline function hydraulic_conductivity(psi, model::VanGenuchten)
    if psi >= zero(psi)
        return model.saturated_hydraulic_conductivity
    end

    S_e = effective_saturation(psi, model)
    return model.saturated_hydraulic_conductivity * S_e^model.pore_connectivity *
           (one(S_e) - (one(S_e) - S_e^inv(model.m))^model.m)^2
end
end # @muladd
