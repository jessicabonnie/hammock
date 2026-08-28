# Mixed-stride v2: register-equality versus inclusion-exclusion

This experiment compares the two hammock similarity outputs on the 20-file
Maurano DHS corpus after `mixed-stride` became the public v2 default.

The register-equality arm is newly timed at p=18, eight threads, five
replicates, and all nine previously studied subB rates:

`1, 0.3, 0.1, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001`.

The inclusion-exclusion (+IE) arm and BEDTools baseline are reused from
`experiments/mixed_stride_v2_figure3_reproduction/`; they were collected on
the same corpus with the same p, thread count, and replicate count. Accuracy
for both similarity types is MAE against exact BEDTools Jaccard over the same
190 unique tissue pairs.

The primary graph selects subB `1`, `0.1`, `0.01`, and `0.001`, with two bars
per rate (register-equality and inclusion-exclusion). The BEDTools median is a
dashed timing reference rather than a third bar.

## Files

- `sbatch_register_equality.sh`: runs the nine-rate register-equality sweep.
- `prepare_inputs.py`: validates and summarizes the raw sweep, then combines
  the selected cells with the existing +IE summary.
- `plot_dual_similarity.R`: renders the paired-bar graph.
- `sbatch_finalize.sh`: prepares data and renders after the sweep succeeds.
- `submit.sh`: submits the benchmark and dependent finalizer.

Generated data, logs, and figures stay inside this directory.
