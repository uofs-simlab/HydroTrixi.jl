@doc raw"""
    ControllerTwoThreshold(semi, indicator;
                           base_level, max_level,
                           coarsen_threshold, refine_threshold)

Create a two-threshold refinement/coarsening AMR controller. Let ``l_k`` be the current
refinement level, let ``l_{\min}`` and ``l_{\max}`` be `base_level` and `max_level`, and
let ``\eta_{\mathrm{c}}`` and ``\eta_{\mathrm{r}}`` be `coarsen_threshold` and
`refine_threshold`, respectively. Elements with ``\eta_k \le \eta_{\mathrm{c}}`` are
coarsened toward ``l_{\min}``, elements with ``\eta_k > \eta_{\mathrm{r}}`` are refined
toward ``l_{\max}``, and all other elements retain their current level.
"""
function ControllerTwoThreshold(semi, indicator;
                                base_level, max_level,
                                coarsen_threshold, refine_threshold)
    return Trixi.ControllerThreeLevel(semi, indicator;
                                      base_level,
                                      med_level = -1,
                                      med_threshold = coarsen_threshold,
                                      max_level,
                                      max_threshold = refine_threshold)
end
