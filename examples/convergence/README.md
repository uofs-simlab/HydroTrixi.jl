# Richards equation convergence studies

Spatial and temporal convergence studies for the pressure-head and mixed forms of the
Richards equation. See the paper for method details.

## Environment and commands

Run from the repository root. Create the local, Git-ignored `run/` environment once:

```sh
mkdir -p run
julia --project=run -e 'using Pkg; Pkg.develop(path="."); Pkg.add(["Trixi", "SciMLBase", "CairoMakie", "LaTeXStrings"])'
```

Run all six studies (60 solves and six PDF figures):

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/convergence/run_richards_convergence.jl
```

For one study, append boundary pair (`DD`/`DN`), refinement (`space`/`time`), and degree:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/convergence/run_richards_convergence.jl DN space 4
```

Replot saved tables without solving; additional directories may be supplied, with each
study present in only one input directory:

```sh
julia --project=run examples/convergence/plot_richards_convergence.jl \
  plots/richards_convergence/RUN_ID
```

For repeated runs, reuse a REPL:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run -i
```

```julia
using LinearAlgebra
BLAS.set_num_threads(1)
include("examples/convergence/run_richards_convergence.jl")
RichardsConvergence.run_convergence(; output_directory = "NEW_RUN_DIRECTORY")
```

## Results and configuration

Error tables are saved in `plots/richards_convergence/RUN_ID/data/*.dat`, and PDF figures
in `plots/richards_convergence/RUN_ID/figures/TIMESTAMP/*.pdf`. Replots use the first input
directory. Earlier and partial tables are preserved; retry a failed run in a fresh directory.

| Study | Boundaries | Degree | Elements | Fixed time step (s) |
| --- | --- | --- | --- | --- |
| Spatial | DD, DN | 3, 4 | 16, 32, 64, 128, 256 | 0.1 |
| Temporal | DD, DN | 3 | 4096 | 4, 2, 1, 0.5, 0.25 |

Both formulations run to 120 s.
