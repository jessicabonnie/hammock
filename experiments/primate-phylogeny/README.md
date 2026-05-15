# primate-phylogeny — Mammalian Phylogeny Recovery from Villar 2015

Tests whether minimizer-HLL sketch similarity over H3K4me3 (or H3K27ac)
peak FASTAs in **adult liver** across 20 mammals reconstructs the known
mammalian phylogeny.

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
snakemake --dryrun
snakemake --profile workflow/slurm_profile/
```

Outputs land at `/vast/blangme2/jbonnie/hammock/primate-phylogeny/results/`
(symlinked here as `results/`).

## Switching marks

To run with H3K27ac instead of H3K4me3, change `mark:` in
`config/config.yaml`. Outputs are namespaced by mark so both can coexist.

## Success criteria

| outcome | what it means |
|---|---|
| Recovered NJ tree has primate + rodent + carnivore + cetartiodactyl clades visible; RF distance to TimeTree significantly below random | **Positive** — sketch recovers phylogeny |
| Pairwise distance correlates with divergence time (Spearman ρ > 0.5) but topology imperfect | **Partial** — distance signal works |
| No structure | **Null** — peak-FASTA sketches don't carry phylogenetic signal in liver |
