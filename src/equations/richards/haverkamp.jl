@muladd begin
#! format: noindent

@doc raw"""
    Haverkamp(; saturated_hydraulic_conductivity, a, beta, b, gamma, theta_s, theta_r)

A Haverkamp constitutive model for the Richards equation,
written in the form used by Celia et al. (1990) and Ireson et al. (2023).
For ``\psi < 0``, the water content and hydraulic conductivity are
```math
\vartheta(\psi) = \theta_{\mathrm{r}} +
\frac{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}
{1 + (a_{\mathrm{H}}|\psi|)^{\beta_{\mathrm{H}}}},
\qquad
\kappa(\psi) = \frac{\kappa_{\mathrm{s}}}
{1 + (b_{\mathrm{H}}|\psi|)^{\gamma_{\mathrm{H}}}},
```
and the effective saturation is
```math
S_{\mathrm{e}}(\psi) =
\frac{\vartheta(\psi)-\theta_{\mathrm{r}}}
{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}.
```
The inverse pressure-head scales ``a_{\mathrm{H}}`` and ``b_{\mathrm{H}}``
correspond to the keywords `a` and `b`.
The water-retention exponent ``\beta_{\mathrm{H}}`` and hydraulic-conductivity exponent
``\gamma_{\mathrm{H}}`` correspond to `beta` and `gamma`.
The parameters ``\kappa_{\mathrm{s}}``, ``\theta_{\mathrm{s}}``, and
``\theta_{\mathrm{r}}`` correspond to `saturated_hydraulic_conductivity`, `theta_s`, and
`theta_r`, respectively.
For ``\psi\geq 0``, the model uses the saturated values
``\vartheta(\psi)=\theta_{\mathrm{s}}`` and ``\kappa(\psi)=\kappa_{\mathrm{s}}``.

# References
- Haverkamp, R., Vauclin, M., Touma, J., Wierenga, P. J., Vachaud, G. (1977).
  A comparison of numerical simulation models for one-dimensional infiltration.
  *Soil Science Society of America Journal*, 41(2), 285-294.
  [DOI: 10.2136/sssaj1977.03615995004100020024x](https://doi.org/10.2136/sssaj1977.03615995004100020024x)
- Celia, M. A., Bouloutas, E. T., Zarba, R. L. (1990).
  A general mass-conservative numerical solution for the unsaturated flow equation.
  *Water Resources Research*, 26(7), 1483-1496.
  [DOI: 10.1029/WR026i007p01483](https://doi.org/10.1029/WR026i007p01483)
- Ireson, A. M., Spiteri, R. J., Clark, M. P., Mathias, S. A. (2023).
  A simple, efficient, mass-conservative approach to solving
  Richards' equation (openRE, v1.0). *Geoscientific Model Development*, 16, 659-677.
  [DOI: 10.5194/gmd-16-659-2023](https://doi.org/10.5194/gmd-16-659-2023)
"""
struct Haverkamp{RealT}
    saturated_hydraulic_conductivity::RealT
    a::RealT
    beta::RealT
    b::RealT
    gamma::RealT
    theta_s::RealT
    theta_r::RealT
end

function Haverkamp(; saturated_hydraulic_conductivity, a, beta, b, gamma, theta_s,
                   theta_r)
    RealT = promote_type(typeof(saturated_hydraulic_conductivity), typeof(a),
                         typeof(beta), typeof(b), typeof(gamma),
                         typeof(theta_s), typeof(theta_r))
    return Haverkamp{RealT}(saturated_hydraulic_conductivity, a, beta, b, gamma,
                            theta_s, theta_r)
end

@inline function effective_saturation(psi, model::Haverkamp)
    if psi >= zero(psi)
        return one(psi)
    end

    return inv(one(psi) + (model.a * abs(psi))^model.beta)
end

@inline function hydraulic_conductivity(psi, model::Haverkamp)
    if psi >= zero(psi)
        return model.saturated_hydraulic_conductivity
    end

    return model.saturated_hydraulic_conductivity /
           (one(psi) + (model.b * abs(psi))^model.gamma)
end
end # @muladd
