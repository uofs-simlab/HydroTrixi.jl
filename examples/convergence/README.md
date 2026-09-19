# Richards equation convergence studies

These examples compare spatial and temporal convergence of the pressure-head
and mixed forms of Richards' equation using a manufactured solution.

## Reproducibility instructions

Run all commands from the repository root. First, create a local `run/`
environment with HydroTrixi and the packages needed for these examples:

```sh
mkdir -p run
julia --project=run -e 'using Pkg; Pkg.develop(path="."); Pkg.add(["Trixi", "SciMLBase", "CairoMakie", "LaTeXStrings"])'
```

The `run/` directory is ignored by Git. Use `--project=run` for the commands below.

Run the full suite: six studies, 60 solves and six PDF figures.

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=run \
  examples/convergence/run_richards_convergence.jl
```

For one study, pass the boundary pair (`DD`/`DN`), refinement (`space`/`time`)
and polynomial degree $N$:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=run \
  examples/convergence/run_richards_convergence.jl DN space 4
```

To redraw figures from saved results:

```sh
julia --project=run examples/convergence/plot_richards_convergence.jl \
  plots/richards_convergence/RUN_ID
```

Pass additional run directories to combine saved studies without solving again.

## REPL usage

For repeated runs, keep the Julia REPL open to reuse package loading and compilation:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=run -i
```

```julia
include("examples/convergence/run_richards_convergence.jl")

# All six studies
RichardsConvergence.run_convergence()

# One study
RichardsConvergence.run_convergence([(bc = "DN", kind = "space", N = 4)])
```

To redraw figures from saved results in the REPL, load the plotting script:

```julia
include("examples/convergence/plot_richards_convergence.jl")

# Plot one saved run
RichardsConvergencePlots.plot_convergence("plots/richards_convergence/RUN_ID")

# Combine studies from multiple saved runs
RichardsConvergencePlots.plot_convergence(
    "plots/richards_convergence/RUN_A",
    "plots/richards_convergence/RUN_B",
)
```

Replace the placeholder run IDs with existing results directories. These calls
use saved tables without solving again. When combining runs, each study must
appear in only one input directory. New PDFs go in a timestamped `figures/`
subdirectory of the first directory passed (`RUN_A` in the second example).

## Saved results

Results go to a new `plots/richards_convergence/RUN_ID/` directory:

| Path within the run | Contents |
| --- | --- |
| `data/` | Error tables |
| `figures/TIMESTAMP/` | PDF figures |

Replots go in the first run directory passed to the plot script. Earlier results are preserved.
The runner prints one line per study and shows warnings. A failure stops the
run and preserves partial tables; retry with a fresh run.

Results are local files, ignored by Git.

## Numerical settings

| Study | Boundaries | Degree $N$ | Elements $K$ | Fixed $\Delta t$ ($\mathrm{s}$) |
| --- | --- | --- | --- | --- |
| Spatial | DD, DN | 3, 4 | 16, 32, 64, 128, 256 | 0.1 |
| Temporal | DD, DN | 3 | 4096 | 4, 2, 1, 0.5, 0.25 |

DD uses penalty Dirichlet conditions at $z=0$ and $z=L$. DN uses penalty
Dirichlet at $z=0$ and the manufactured Neumann flux
$f(\psi,\partial_z\psi)$ at $z=L$. The penalty factor is 1 in both cases.

Every study runs both formulations on a column of length $L=0.2\,\mathrm{m}$
until $t=120\,\mathrm{s}$, using Float64, DGSEM/LDG and Rodas5P/KLU. Sparse
forward AD supplies a fresh Jacobian each step. Time steps are fixed; explicit
time stops prevent extra steps from roundoff. The mesh is not adapted during a
solve.

For a uniform mesh with $K=2^l$ elements, the nominal mesh size is

$$
\Delta z = \frac{L}{KN}.
$$

Spatial refinement uses levels $l=4,\ldots,8$. Errors in pressure head $\psi$
are measured in metres: the length-normalized $L^2$ error uses collocated LGL
quadrature with $N+1$ nodes per element, and the nodal $L^\infty$ error is the
maximum absolute error over those nodes. The tables
also contain DOFs, observed convergence orders, solver status and accepted step
counts.

Figures show all five refinement levels and both error norms. Reference
triangles mark order $N+1$ for spatial refinement and orders 4 and 5 for
temporal refinement; measured orders are listed in the tables.
