# Corrected cross-reference Mode D rerun at p=18

This is the precision-18 companion to `../rehash_rerun/`. It repeats the exact
63-cell parameter union, both peak types, and seeds 1, 42, and 99 using the
archived FASTAs and corrected `rehash-selector64` sequence-HLL path. It changes
only HLL precision from 24 to 18.

Raw matrices and Slurm logs are isolated under:

```
/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64_p18/
```

Submit from the repository root:

```bash
mkdir -p /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64_p18/logs/slurm
sbatch --array=0-377%50 experiments/ref-comparison/rehash_rerun_p18/run_array_task.sh
python experiments/ref-comparison/rehash_rerun/analyze.py \
  --results /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64_p18 \
  --cells experiments/ref-comparison/rehash_rerun/cells.tsv \
  --metadata docs/data/exp_a_metadata.tsv \
  --output-dir /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64_p18/summaries \
  --precision 18
```

The manuscript's clade claim at $k\geq8$ is checked directly by clustering every
matrix (`results/clade_recovery.csv`; requires SciPy). At $k=5$, results flagged in
`tied_merge_heights` can depend on the clustering implementation.

```bash
python experiments/ref-comparison/rehash_rerun/verify_clades.py \
  --results /vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64_p18 \
  --cells experiments/ref-comparison/rehash_rerun/cells.tsv \
  --output experiments/ref-comparison/rehash_rerun_p18/results/clade_recovery.csv \
  --precision 18
```
