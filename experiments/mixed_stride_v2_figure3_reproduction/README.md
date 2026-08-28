# Figure 3 reproduction with mixed-stride v2 as default

This directory refreshes Figure 3B in `paper/draft.md` without overwriting
canonical manuscript data or `paper/figures/pairwise_scaling.png`. Figure 3A
is unchanged and is rendered from its existing canonical results.

## Protocol

Panel A reuses `docs/data/cpp_vs_bedtools_t16_p18.csv`, the existing seeded
synthetic scaling benchmark at p=18 and 16 threads:

- 10,000 intervals per BED file;
- N = 2, 4, 8, 16, 32 with 20 replicates;
- N = 64, 128, 256, 512 with 3 replicates;
- corpus seed 20260810;
- bedtools 2.30.0, hammock's register-equality arm, and the recommended
  `jaccard_similarity_ie` arm used in the main-text figure.

Panel B repeats the 20-file Maurano benchmark at p=18 and eight threads:

- public/default `mixed-stride` at subB = 1.0, 0.3, 0.1, 0.03, 0.01,
  0.008, 0.005, 0.003, and 0.001;
- five hammock replicates per rate, emitting `jaccard_similarity_ie`;
- five bedtools 2.30.0 batches over the 190 unique pairs;
- MAE computed against exact bedtools Jaccard.

At subB = 1.0, 0.1, 0.01, 0.008, 0.005, and 0.001 the reciprocal stride is
integral, so this branch's v2 default delegates to the legacy grid. The 0.3,
0.03, and 0.003 cells exercise the new within-chromosome mixed-gap behavior.

## Files

- `sbatch_panel_b.sh`: Maurano subsampling and bedtools job.
- `prepare_inputs.py`: summarizes Panel B for the plotting script.
- `sbatch_finalize.sh`: prepares inputs and renders the reproduced figure.
- `submit.sh`: submits the Panel B benchmark and a dependent finalizer.

All generated files land under `results/`, `logs/`, and `figures/` here. The
finalizer renders `pairwise_scaling_reproduced.png` with the paper's original
three Panel B rates and `pairwise_scaling_reproduced_expanded.png` with all
nine rates. Panel A is identical in both figures and does not use subB.
