# AGENTS.md

## Formatting
- Use the project `JuliaFormatter.jl` configuration, including the 92-column margin.
- Wrap docstrings and comments as close to 92 characters as practical. Leave equations,
  DOI links, and similar constructs unwrapped when tighter wrapping would reduce clarity.
- Do not introduce new identifiers with leading underscores in tracked code.
- Do not preserve backward compatibility at all.
- Follow Trixi.jl style and conventions.
- In prose, always write "Trixi.jl" and "HydroTrixi.jl" with the ".jl" suffix.
- Write "pressure head" without a hyphen when referring to the quantity. Use
  "pressure-head" only as a compound adjective, such as "pressure-head form".
- Unless at the very end, run the example elixirs, not the full test suite.
- Use `julia --project=run` for plotting examples and scripts with plotting dependencies.

## Dependencies
- Adding dependencies is acceptable only after obtaining explicit approval from the user.

## Validation
- Add explicit validation only when invalid input could otherwise run successfully and
  silently produce incorrect results. Let unsupported signatures and inputs that will
  inevitably fail later use ordinary dispatch, bounds, property, or similar errors instead
  of adding pre-emptive checks solely to provide a friendlier exception.
