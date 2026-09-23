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

This runs 18 simulations and creates six PDF figures. To run one three-tolerance
study, pass its benchmark and case:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/conservation/run_richards_conservation.jl haverkamp mixed
```

Benchmarks are `haverkamp` and `new_mexico`. Cases are `mixed`,
`pressure_head`, and `pressure_head_water_content_transfer`.

To redraw figures without solving again, pass one or more saved run directories:

```sh
julia --project=run examples/conservation/plot_richards_conservation.jl \
  plots/richards_conservation/RUN_ID
```

Results are written to a new `plots/richards_conservation/RUN_ID/` directory.
The `data/` directory contains accepted-step histories, while timestamped
subdirectories of `figures/` contain the PDFs. Earlier and partial results are
preserved.

## Numerical settings

All cases use polynomial degree 3, initial mesh level 6, sparse forward AD, and
the default Rodas5P/KLU algorithm. The normalized water-content TV indicator is
applied every 10 accepted steps with base level 2, maximum level 10, and
coarsen/refine thresholds 0.003/0.03. Initial AMR only refines.

Time stepping starts from 0.01 s without a positive minimum-step floor and
controls water-content error using `abstol = 1e-11` and
`reltol = 1e-5, 1e-7, 1e-9`. The mixed form transfers water content under AMR.
The pressure-head form is run once transferring pressure head and once
transferring water content.

Mass bias is recorded after every accepted step. Figures show its absolute
magnitude on common logarithmic limits; the saved tables retain its sign.
The plotting script reads only these tables and can regenerate every figure
without constructing or solving the semidiscretizations.
