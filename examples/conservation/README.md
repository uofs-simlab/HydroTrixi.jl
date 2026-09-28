# Richards equation conservation studies

These examples compare fully discrete water mass bias for the Celia-Haverkamp
and Celia-New Mexico infiltration problems with adaptive meshes and time steps.

## Running the studies

Run from the repository root using the `run/` environment described in
`examples/convergence/README.md`:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/conservation/run_richards_conservation.jl
```

This runs 24 simulations and creates four grouped mass-bias PDFs, four grouped
time-step-history PDFs, and 48 solution-mesh PDFs (half-time and full-time for
each solve). Water-content and pressure-head transfer are overlaid in the grouped
figures. To run one three-tolerance study, pass its benchmark and case:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/conservation/run_richards_conservation.jl haverkamp mixed
```

Benchmarks are `haverkamp` and `new_mexico`. Cases are `mixed`,
`mixed_pressure_head_transfer`, `pressure_head`, and
`pressure_head_water_content_transfer`.

To redraw figures without solving again, pass one or more saved run directories:

```sh
julia --project=run examples/conservation/plot_richards_conservation.jl \
  plots/richards_conservation/RUN_ID
```

Results are written to a new `plots/richards_conservation/RUN_ID/` directory.
The `data/` directory contains accepted-step histories. The `snapshots/`
directory contains plot-ready numerical pressure-head and mesh-edge tables at
exactly half-time and full-time. Timestamped subdirectories of `figures/`
contain the PDFs. Replotting reads these tables and creates a new figure
directory. Older runs without snapshots still produce mass-bias and time-step
PDFs; the script reports which snapshot PDFs it skipped. Earlier and partial
results are preserved.

## Numerical settings

All cases use polynomial degree 3, initial mesh level 6, sparse forward AD, and
the default Rodas5P/KLU algorithm. The normalized water-content TV indicator is
applied every 10 accepted steps with base level 2, maximum level 10, and
coarsen/refine thresholds 0.003/0.03. Mesh adaptation starts after the initial state.

Time stepping starts from 0.01 s without a positive minimum-step floor and
controls water-content error using `abstol = 1e-11` and
`reltol = 1e-5, 1e-7, 1e-9`. Each formulation is run once transferring water
content and once transferring pressure head under AMR.

Mass bias is recorded after every accepted step. Figures show its absolute
magnitude on common logarithmic limits; the saved tables retain its sign.
Time-step figures plot each accepted step's size against time on a linear
axis and omit the initial step-size guess. Solution-mesh PDFs show pressure
head against depth with vertical element-boundary guides, using the same figure
style as the AMR animations.
Colors distinguish tolerances, while solid and dashed lines indicate water-content
and pressure-head transfer, respectively. Grouped figures have a three-entry
tolerance legend. A single-case run plots its available transfer choice.
The plotting script can regenerate every available figure without constructing
or solving the semidiscretizations.
