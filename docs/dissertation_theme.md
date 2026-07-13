# Dissertation Theme Discussion (Working Notes)

## Purpose

These notes are a living design document for developing the overarching dissertation narrative. They are intentionally exploratory rather than polished prose. The goal is to identify a central scientific question that naturally connects the research chapters instead of forcing three independent papers into a single dissertation.

## Current Status

The dissertation currently consists of three research chapters:

- Chapter 2: DandD (published)
- Chapter 3: Hammock (currently being written)
- Chapter 4: RNA-seq genotype inference for sequence-to-function models (planned)

The objective is to develop a coherent scientific story spanning all three chapters.

---

## Candidate Theme A

### Computational methods for extending the utility of large-scale genomic resources

Central question:

> How can computational methods allow existing genomic resources to continue producing new biological insight as genomics evolves?

The motivation is that genomic resources continue to accumulate while sequencing technologies, reference genomes, and computational methods evolve.

---

## Candidate Theme B (currently preferred)

### Representational methods for extending the utility of genomic resources

Each chapter develops a new representation of existing biological information that enables analyses that were previously impractical or impossible.

Representation is viewed as the mechanism. Extending the scientific utility of genomic resources is the objective.

---

## Emerging refinement

A stronger framing may be that the dissertation is fundamentally about **adapting biological information to new analytical contexts**.

Across genomics, biological measurements often outlive the technologies, reference genomes, and computational methods that were available when they were generated.

Rather than replacing these datasets, the dissertation develops computational approaches that adapt existing biological information so it remains useful for answering new scientific questions.

Possible central question:

> How can biological information be adapted to remain scientifically useful as analytical methods, reference genomes, and computational capabilities evolve?

Representations become one mechanism of adaptation rather than the ultimate objective.

---

## Research philosophy

None of the chapters require collecting new biological data.

Instead, each develops computational methods that allow existing biological information to answer questions that were previously difficult or impossible.

The emphasis is therefore on increasing the scientific value extracted from existing data rather than generating new data.

---

## Longer-term research direction

One recurring personal research interest is that carefully chosen **representations can reveal fundamental characteristics of biological datasets** while simultaneously enabling new analyses.

Examples:

- δ summaries expose redundancy and diversity in expanding genome collections.
- Interval sketches preserve biologically meaningful similarity while enabling scalable comparison.
- Inferred genotypes recover latent genomic variation that enables modern sequence-aware models.

This may represent a broader research philosophy beyond the scope of this dissertation. It currently fits Chapters 2 and 3 more naturally than Chapter 4, so it should be treated as an emerging idea rather than the primary dissertation theme.

---

## Working observation

A recurring pattern across all chapters is:

Existing biological resource → New scientific question → Original representation becomes insufficient → Computational adaptation → Resource becomes useful in a new analytical context.

This pattern may ultimately provide the strongest unifying narrative for the dissertation.

---

## Open questions

We have not yet decided whether the dissertation should primarily emphasize:

1. Computational methods.
2. Representational methods.
3. Adapting biological information to new analytical contexts.
4. Extending the utility of genomic resources.
5. Some combination of these.

Another unresolved question is whether revealing fundamental characteristics of biological datasets through representation is a dissertation theme or a broader long-term research philosophy.