@doc raw"""
    ControllerTwoThreshold(semi, indicator;
                           base_level, max_level,
                           coarsen_threshold, refine_threshold)

Create a two-threshold refinement/coarsening AMR controller.
Let ``l_k`` be the current refinement level,
let ``l_{\min}`` and ``l_{\max}`` be `base_level` and `max_level`,
and let ``\eta_{\mathrm{c}}`` and ``\eta_{\mathrm{r}}``
be `coarsen_threshold` and `refine_threshold`, respectively.
Mark an element for refinement when ``l_k<l_{\max}`` and ``\eta_k(t)>\eta_{\mathrm{r}}``,
and for coarsening when ``l_k>l_{\min}`` and ``\eta_k(t)\leq\eta_{\mathrm{c}}``.
All other elements retain their current level.
The AMR callback merges a pair of siblings only when both are marked for coarsening;
the tree adaptation also enforces neighbouring refinement levels that differ by at most one.
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
