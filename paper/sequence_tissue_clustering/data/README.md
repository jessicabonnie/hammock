Paper-local inputs for the Maurano sequence-mode tissue figures. The current
single-run and precision-comparison candidates use prespecified seed 42; the
paired-cut diagnostic uses representative seed 12. Each active CSV below is a
byte-for-byte copy of a validated independent Mode D
`rehash-selector64` output produced from source commit
`008e9f24653cf4f40d10305c3b56ef29979ea6ae`:

| Local input | SHA-256 |
|---|---|
| `p12_seed00042_k10_w30.csv` | `5c7014c64295b0dec60387e7f4d78f9d0a9b9e390a1406d7a7c59800fcc50e73` |
| `p18_seed00042_k10_w30.csv` | `de146ace640cfa153e1dc63c4bcb17e48f10597f0b641951d193e14fd6a9e179` |
| `p24_seed00042_k10_w30.csv` | `bd0b42953d0bd44207f6b3f5668dba6cfabd2eca75f24c13193084a134b6ffe2` |
| `p18_seed00042_k20_w20.csv` | `1ff426bf2a092b4cd14de24ccc1bc53f729b111ccea9c5962e3f903e6bcf3434` |
| `p18_seed00012_k20_w20.csv` | `6489cceb0e207f8c6d651b6fd358bdefb2103422afbfd6eb86efc47841f854bc` |

The older seed-0 and seed-1 paper-local copies are retained for reproducing
earlier drafts, but are no longer the defaults.

Run `paper/sequence_tissue_clustering/render_rehashed_tissue_figures.sh` from
the repository root after loading R 4.3.0. It regenerates the ten- and
eight-class single-run dendrogram/heatmap pairs, both fixed-seed precision
composites, and the seed-12 paired-cut supplementary diagnostic. The render
path writes vector PDFs first and rasterizes them at 300 dpi with `pdftoppm`,
avoiding a dependency on the unavailable R Cairo binary.
Run the same command with `--verify-only` to check all five pinned rehash-input
hashes without rendering.

`p18_seed00012_k20_w20.csv` supports the paired eight- and ten-cluster-cut
dendrogram. Seed 12 was selected after the 101-seed extension as the run with
the smallest mean absolute deviation from the elementwise across-seed median
distance matrix over the 190 off-diagonal sample pairs (MAE 0.000887657). The
same matrix is shown in both panels. Its lowest average-linkage merge is at
1-J = 0.367149, so the display uses an axis break at 0.35; the break therefore
compresses terminal stem space without compressing any observed merge. The
display rotates equivalent dendrogram branches by annotated tissue rank so the
two leg samples are adjacent; assertions verify that cophenetic distances and
merge heights are unchanged.
