# Corrected cross-reference Mode D rerun

This rerun applies `--sequence-hll-hash rehash-selector64` to the exact 63-cell
parameter union behind Figure 5. It reuses the archived broad- and narrow-peak
FASTAs and does not repeat alignment, peak calling, or sequence extraction.

The primary paper-facing realization is seed 42. Seeds 1 and 99 are frozen
sensitivity replicates. All runs use Mode D, precision 24, one sketching thread,
and the full metric output. The 378 array tasks are ordered by seed, peak type,
then the rows of `cells.tsv`.

Results are written separately under:

```
/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/
```

Submit and analyze from the repository root:

```bash
mkdir -p /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/logs/slurm
sbatch --array=0-377%50 experiments/ref-comparison/rehash_rerun/run_array_task.sh
python experiments/ref-comparison/rehash_rerun/analyze.py \
  --results /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64 \
  --cells experiments/ref-comparison/rehash_rerun/cells.tsv \
  --metadata docs/data/exp_a_metadata.tsv \
  --output-dir /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/summaries
```
