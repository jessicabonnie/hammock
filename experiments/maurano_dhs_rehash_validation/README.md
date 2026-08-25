# Maurano Mode D independent-rehash validation

This experiment asks whether the historical ten-class Mode D optimum at
`k=10,w=30` survives exact selected-feature comparison and independent 64-bit
HLL rehashing. It deliberately separates three targets: exact selector-feature
Jaccard (estimator truth), ten-class tissue recovery (an exploratory biological
objective), and BEDTools coordinate Jaccard (an external application target).
The prespecified ten-class endpoint remains primary. A secondary eight-class sensitivity endpoint
collapses `fMuscle_arm`, `fMuscle_back`, and `fMuscle_leg` into `fMuscle` and cuts the same average-linkage
hierarchy into eight groups. It is scored seed-for-seed from the completed matrices; it does not add HLL runs.

## Frozen design

The primary grid contains the 37 valid cells from `k={8,10,15,20,25}` and
`w={8,10,20,30,50,100,200,300,500}`, retaining `w>=k`. Rehashed HLL runs use
`p=18`, seeds `1,2,3,17,42,43,99,31337`, `--threads 1`, and the production
`rehash-selector64` implementation. The exact experiment uses the same
`digest.window_minimizer` selector path and represents selector identities as
uint32 values; no-minimizer whole-record fallbacks are kept as literal,
upper-cased byte strings in a separate feature domain.

The exact ranking sorts by ten-class ARI, then NMI, and reports unresolved ties.
The rehash ranking uses median tissue ARI over all prespecified seeds, full ARI
range, exact-partition recovery frequency, median NMI, and median MAE against
exact Jaccard. No result is selected from the best observed seed. Precision
follow-up at `p={12,18,22,24}` is generated only after p=18 leaders are frozen,
for historical `k=10,w=30`, every exact leader, and any distinct robustness
leader.

An exploratory extension evaluates 101 total seeds (0 through 99 plus 31337)
at the historical `k=10,w=30` cell for all four precisions. The original eight
seeds remain the frozen confirmatory phases; only the other 93 seeds are run by
`--phase extension`. Exact-partition frequencies include Wilson 95% intervals.
An additional interpolation phase evaluates the missing p=13–17, p=19–21,
and p=23 values over all 101 seeds, completing the p=12–24 series. It reports
seed-to-seed ARI/NMI dispersion and discrete partition and hierarchy-state
frequencies, rather than relying on medians alone.

Because parameter selection and recovery use the same 20 labelled samples,
all tissue rankings are exploratory and are not estimates of generalization.

## Layout and commands

All commands are run from this directory with the isolated build's Python:

```bash
python scripts/inventory_inputs.py --config config.yaml
python scripts/extract_exact_features.py --config config.yaml --dry-run
python scripts/run_rehash_sweep.py --config config.yaml --hammock /path/to/isolated/hammock --dry-run
python scripts/validate_outputs.py --config config.yaml --all
python scripts/analyze_rehash_sweep.py --config config.yaml
```

`inventory_inputs.py` refuses a grid other than the historical 37 cells and
writes checksummed input/provenance manifests. Exact extraction is resumable at
the sample/cell level. The HLL runner stages each command in a temporary
directory, validates the emitted CSV by column name, and atomically installs it
with a completion manifest. Existing partial or invalid outputs are moved to a
quarantine directory rather than overwritten. Bounded Slurm entry points live
under `workflow/`, but submission remains a separate, explicitly authorized
step after review of the dry-run job table.

## Correctness gates

Every matrix must have 20 samples, 400 ordered pairs, finite bounded values,
symmetry, and a unit diagonal. Repeated exact extraction must be byte-identical.
Fixed configuration/seed HLL reruns must be deterministic after timing and path
provenance are normalized. Seed and hash-mode metadata are included in every
output name and completion manifest. Exact fixtures compare production selector
extraction with direct Python set operations.

The full biological sweep is substantial: 37 exact extractions over 20 large
FASTAs plus 296 p=18 HLL jobs, followed by a leader-dependent precision sweep.
Do not launch it on the scheduler without explicit approval. A two-sample smoke
test and one complete 20-sample cell must pass first, and measured peak RSS must
be used to set concurrency.

## Interpretation

The final decision must classify `k=10,w=30` as exactly one of: exact and
robust; exact but estimator-sensitive; part of an exact plateau; not the
exact-feature optimum; or unresolved. A single seed is never sufficient to
change the compatibility default or manuscript conclusion.
