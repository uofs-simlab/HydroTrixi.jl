# Repository Agent Notes

## Tests

A systematic unit-test framework will be introduced later. For now, add coverage
through the existing Trixi.jl-style elixir regression tests in `test/runtests.jl`,
using `@trixi_testset` and `@test_trixi_include` with reference `l2` and `linf` errors
where appropriate. Do not add ad hoc unit checks or standalone unit-test files.
Preserve existing tests unless a change is explicitly requested.

## Convergence results

Save numerical tables and figures only. Do not generate metadata, provenance,
source/environment snapshots, plot revision logs, question-based reports,
Markdown run reports, experiment-index entries, or timing measurements.
Preserve earlier results and partial tables from failed runs. Follow the
convergence README; do not silently change numerical settings after a failure.

## Reusing a Julia session

For repeated simulations or animations, keep a Julia REPL open so Julia startup,
package loading, and compilation are reused between runs.

1. From the repository root, start Julia in a PTY and retain the returned session ID:

   ```sh
   julia --project=run -i
   ```

   Set up the local `run/` environment as described in
   `examples/convergence/README.md`; it includes CairoMakie and LaTeXStrings for
   visualization. For non-visual core simulations, `julia --project=. -i` is sufficient.

2. Send commands to that process with `write_stdin`. For example:

   ```julia
   redirect_stdout(devnull) do
       include("examples/visualization/animate_richards_celia_new_mexico_amr.jl")
   end
   println("RUN_COMPLETE")
   ```

   Redirecting standard output avoids retaining the verbose time-step log; exceptions
   and other standard-error output remain visible. Omit the redirection when detailed
   solver diagnostics are needed.

3. Poll the same session with an empty `write_stdin` call until the Julia prompt returns.
   To rerun after editing a file, send the same `include(...)` command to that session.
   Do not send `exit()` or Ctrl-D until all requested iterations are complete.

Julia processes in this repository may need to run outside the filesystem sandbox because
MPI initialization opens a local listener. Request the corresponding execution approval
when starting the persistent REPL.
