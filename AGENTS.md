# Repository Agent Notes

## Tests

A systematic unit-test framework will be introduced later. For now, add coverage
through the existing Trixi.jl-style elixir regression tests in `test/runtests.jl`,
using `@trixi_testset` and `@test_trixi_include` with reference `l2` and `linf` errors
where appropriate. Do not add ad hoc unit checks or standalone unit-test files.
Preserve existing tests unless a change is explicitly requested.

Independent tests and diagnostic simulations may run concurrently in multiple Julia
processes, up to six at a time, with Julia and OpenBLAS threads both set to one.
Follow the session reuse and concurrency guidelines below.

## Code style

Use explicit `if` statements for control flow. Do not use short-circuit expressions to
trigger exceptions, returns, loop control, assignments, or side effects, such as
`condition || throw(...)`, `condition && return`, or `condition || mkpath(...)`.
Short-circuit operators remain appropriate within ordinary boolean expressions.

Use ordinary multiple dispatch when changed argument types already select the intended
method. Use `invoke` only when deliberately bypassing a more-specific applicable method;
do not use it merely to delegate after unwrapping an argument to its base type.

## Convergence results

Save numerical tables and figures only. Do not generate metadata, provenance,
source/environment snapshots, plot revision logs, question-based reports,
Markdown run reports, experiment-index entries, or timing measurements.
Preserve earlier results and partial tables from failed runs. Follow the
convergence README; do not silently change numerical settings after a failure.

## Reusing Julia sessions and running concurrently

For repeated tests, simulations, or animations, keep Julia REPLs open so Julia startup,
package loading, and compilation are reused between runs. Reuse an idle compatible
session before starting another one.

1. From the repository root, start Julia in a PTY and retain the returned session ID:

   ```sh
   JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run -i
   ```

   Set up the local `run/` environment as described in
   `examples/convergence/README.md`; it includes CairoMakie and LaTeXStrings for
   visualization. For non-visual core simulations, use `--project=.` instead. In each
   session, explicitly set OpenBLAS threads after loading LinearAlgebra:

   ```julia
   using LinearAlgebra
   BLAS.set_num_threads(1)
   ```

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
   To rerun after editing an elixir or diagnostic script, send the same `include(...)`
   command to that session. Re-including an elixir does not reload changes to package
   source already loaded by `using HydroTrixi`; restart affected sessions or use an
   established source-reloading workflow before validating package changes.
   Keep reusable sessions open until the requested iterations are complete, unless a
   restart is needed for source changes or to release resources.

### Concurrent session guidelines

- Run up to six independent Julia workloads at once when useful. This is a ceiling,
  not a target; use fewer when memory pressure or compilation contention warrants it.
  Count test subprocesses and workers toward the limit, not just interactive REPLs.
- Set `JULIA_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` in the launch environment and
  use `--threads=1`. Keep these settings for test subprocesses as well; do not override
  them with multi-threaded test launch options.
- Keep track of each session ID, active project, loaded source version, and assigned
  workload in the working conversation. Give each session one command stream; never
  send another run into a busy session or let multiple agents drive the same REPL.
- Reuse each session for compatible follow-up runs. Reconstruct mutable problem,
  mesh, integrator, and callback state for independent cases so previous runs do not
  contaminate later ones.
- Give concurrent runs distinct output paths, or disable file output. Do not run
  copies of tests that write the same files concurrently. Preserve existing results.
- Set up dependencies and finish initial precompilation before launching parallel
  workloads. Do not concurrently modify the same package environment with `Pkg`.
- `Pkg.test()` launches a separate test process; keeping its parent REPL open does
  not preserve that child's compilation. For repeated targeted runs, reuse a REPL
  with the required test dependencies and the existing elixir regression framework.

Julia processes in this repository may need to run outside the filesystem sandbox because
MPI initialization opens a local listener. Request the corresponding execution approval
when starting the persistent REPL.
