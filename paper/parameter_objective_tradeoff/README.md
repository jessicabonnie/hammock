The sequence-mode parameter objective tradeoff figure reads the validated
`rehash-selector64` run table at
`experiments/maurano_dhs_rehash_validation/results/summaries/rehash_run_scores.csv`.
It filters to p=18 and the eight prespecified independent seeds (1, 2, 3, 17,
42, 43, 99, and 31337), requires all 37 parameter cells to contain exactly one
run per seed, and plots the per-cell medians of `bedtools_pearson` and either
`ari_10class` or `ari_8class`.

Both the ten- and eight-cut figures encode window size by color and k-mer size
by open or line-only point shapes. In the ten-cut figure, the numerical optimum
is colored magenta and the biological optimum bright orange using their native
shapes. Magenta therefore denotes the numerical optimum in both class views.

The optional fourth argument sets the lower Pearson display limit. It produces
a zoomed view without removing lower-agreement cells from the analysis. For
example, passing `0.9` displays Pearson r=0.90--1.00. Zoomed views omit the
non-plateau points while leaving them in the analysis. Their color limits,
color breaks, and shape levels remain fixed to those of the full graph. Both
eight-cut views color the numerical-optimum point magenta. The full view colors
every tied biological-optimum point bright orange, while the zoom leaves those
points in their normal window-size colors.

Exact selected-feature scores are not inputs to this figure.
