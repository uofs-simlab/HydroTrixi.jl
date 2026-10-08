# Richards equation conservation studies

Water mass-bias studies for the Celia-Haverkamp and Celia-New Mexico problems.
See the paper for method details.

## Run and replot

Run from the repository root using the environment in
[the convergence README](../convergence/README.md). The full suite runs 24 solves:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/conservation/run_richards_conservation.jl
```

For one study at all three tolerances:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --threads=1 --project=run \
  examples/conservation/run_richards_conservation.jl haverkamp mixed
```

Benchmarks: `haverkamp`, `new_mexico`. Cases: `mixed`, `mixed_pressure_head_transfer`,
`pressure_head`, `pressure_head_water_content_transfer`.

Replot saved results without solving; additional directories may be supplied:

```sh
julia --project=run examples/conservation/plot_richards_conservation.jl \
  plots/richards_conservation/RUN_ID
```

## Saved results

Each run saves:

- Accepted-step histories: `plots/richards_conservation/RUN_ID/data/*.dat`.
- Half-time and full-time pressure-head and mesh-edge snapshots:
  `plots/richards_conservation/RUN_ID/snapshots/*.dat`.
- PDF figures: `plots/richards_conservation/RUN_ID/figures/TIMESTAMP/*.pdf`.

Runs use fresh directories; existing and partial results are preserved.
In a persistent REPL, include the runner and call
`RichardsConservation.run_conservation(; output_directory = "NEW_RUN_DIRECTORY")`
to choose a fresh directory.

Replotting creates new figures in the first input directory; missing snapshot tables are
reported and skipped.
Mass-bias plots compare both transfer strategies. Time-step plots use water-content
solution transfer at all three tolerances, as in the manuscript; they are skipped when
only pressure-head-transfer results are supplied.
