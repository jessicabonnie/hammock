# primate-phylogeny — Mammalian Phylogeny Recovery from Villar 2015

Tests whether minimizer-HLL sketch similarity over H3K4me3 (or H3K27ac)
peak FASTAs in **adult liver** across 20 mammals reconstructs the known
mammalian phylogeny.

> **Outcome: PARTIAL, accepted 2026-05-14. Do not read this as a success.**
> Seven of the planned 20 species were run (hsa, rheMac, calJac, Mmus,
> canFam, bosTau, monDom). The primate sub-tree `((hsa, rheMac), calJac)`
> is recovered cleanly and survives every variation tried — (k, w) sweep,
> HLL precision 24/26/28, H3K4me3 ↔ H3K27ac, Jaccard vs cosketch. **The
> deeper topology fails**: Mmus pairs with monDom (long-branch attraction),
> canFam and bosTau do not form a Laurasiatheria clade, and monDom is not
> recovered as the marsupial outgroup. No precision, mark, or metric fixes
> it. The diagnosed bottleneck is the dynamic range of the true Jaccard
> signal (cross-species similarities span ~0.055 at k=10/w=10), not sketch
> noise. The follow-on direction — repeat-masking the peak FASTAs,
> restricting to non-promoter peaks, or k-mer weighting — is unblocked but
> **not done**. Full readout: `docs/experiment_design.md`; pickup brief:
> `memory/project_primate_phylogeny_substrate_engineering.md`.

> **Cross-species Mode D caveat.** Mode D Jaccard between FASTAs from
> different reference genomes measures **shared k-mer content**, which at
> mammalian divergence is dominated by repeats, low-complexity sequence and
> conserved promoter motifs — it is *not* a homology measure
> (`CLAUDE.md`, "Cross-reference caveat"). hammock's default `k=8` is
> unsuitable here; this sweep runs k up to 20 for that reason, and k=5
> cells are reported as saturated/uninformative.

Data: Villar et al. 2015 (*Cell* 160:554), ArrayExpress `E-MTAB-2633`.
Single tissue (liver), 20 species, two histone marks (H3K27ac, H3K4me3)
called by MACS against Ensembl v73 (Sept 2013) reference assemblies.

## Sample design

20 species spanning 7 clades (primates × 4, rodents × 4, lagomorph,
carnivores × 3, cetartiodactyls × 5, scandentia, marsupials × 2).
Full list with peak URLs in `config/samples.tsv`.

## Pipeline shape

```
Villar peak BED (curl from ArrayExpress)
        ↓
bedtools getfasta (per-species v73 reference)
        ↓
hammock all-vs-all (Mode D, (k,w) sweep)
        ↓
neighbor-joining tree + distance matrix
        ↓
(future) topology comparison to TimeTree reference
```

(k, w) sweep: `k ∈ {5, 8, 10, 15, 20}` × `w ∈ {5, 8, 10, 15, 20, 30}` filtered to `w ≥ k` (20 cells).

## Sample inclusion

The Snakefile only runs species whose `ucsc_assembly` key is listed under
`references:` in `config/config.yaml`. This lets us iterate from a small
plumbing run (human + mouse, both have local FASTAs) up to the full
20-species panel as v73-era references are staged.

## Running

```bash
conda activate claude-ref-comparison
ml r/4.3.0                                   # tree building is scripts/build_phylogeny.R
snakemake --dryrun
snakemake --profile workflow/slurm_profile/
# or, as one batch job:
sbatch sbatch_workflow.sh
```

Outputs land at `/vast/blangme2/jbonnie/hammock/primate-phylogeny/results/`
(symlinked here as `results/`), namespaced by mark: `results/{H3K4me3,H3K27ac}/k{k}_w{w}/`.

Other scripts in this directory, none of them wired into the Snakefile:

| script | what it does |
|---|---|
| `scripts/fetch_v73_references.sh` | stages the Ensembl release-73 reference FASTAs the config points at |
| `scripts/precision_probe.sh` | sbatch sweep of HLL precision (24/26/28) at k∈{8,10}, w=10; output under `results/H3K4me3/precision_probe/` |
| `estimator_ie_topology.py` | re-scores the recovered topology on `jaccard_similarity_ie` (derived from the shipped containment columns — no re-sketching). Run as `python estimator_ie_topology.py [--mark H3K4me3 …]`. Its conclusion is not recorded in these docs. |

### Re-running vs. what is on disk

The archived CSVs and trees were produced **before hammock v0.6.0**, so they
carry the `_with_ends` column family and each cell holds six trees
(`{jacc,cosmax,cosgeom}_{mz,we}`). v0.6.0 removed that family
(`CLAUDE.md` divergence #8), so a re-run today produces only the three `_mz`
trees. `config/config.yaml`'s `primary_sim_col` has already been moved to
`jaccard_similarity`. Parse these CSVs **by column name** — the schema
changed and positional field indices do not line up across versions.

## Switching marks

To run with H3K27ac instead of H3K4me3, change `mark:` in
`config/config.yaml`. Outputs are namespaced by mark so both can coexist.

## Success criteria

| outcome | what it means |
|---|---|
| Recovered NJ tree has primate + rodent + carnivore + cetartiodactyl clades visible; RF distance to TimeTree significantly below random | **Positive** — sketch recovers phylogeny |
| Pairwise distance correlates with divergence time (Spearman ρ > 0.5) but topology imperfect | **Partial** — distance signal works |
| No structure | **Null** — peak-FASTA sketches don't carry phylogenetic signal in liver |
