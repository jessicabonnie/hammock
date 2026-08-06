# Working Methods

This document is the working space for determining the content and organization of the Methods section. The current repository README is the source of truth for the implemented software. Decisions made here should be propagated to `paper/outline.md` once agreed.

## Agreed reporting decisions

- The manuscript will not mention the former `with_ends` sequence metrics.
- Inclusion–exclusion Jaccard (`jaccard_similarity_ie`) is the principal interval-mode metric for comparison with `bedtools jaccard`.
- Register equality (`jaccard_similarity`) is retained by the software for compatibility and may be described in Methods, but it should not be the primary validation target in the main text.
- Figure 4 should emphasize inclusion–exclusion Jaccard as the metric that lies near the identity line with BEDTools. Register equality should be removed from the main figure unless it serves a clearly necessary secondary purpose.

# Proposed Methods organization

## 4.1 Software implementation

### 4.1.1 Software architecture and distribution

Include:

- Hammock is distributed as a Python package exposing the `hammock` command-line entry point.
- Python handles command-line orchestration, input-list processing, mode selection, reference handling, BED-to-FASTA workflow management, and tabular output.
- Performance-critical BED parsing and HyperLogLog sketch construction are implemented in a C++17 extension exposed through pybind11.
- The extension is built with `scikit-build-core`, CMake, and a C++17 compiler, with OpenMP support for interval-mode parallelism.
- Hammock uses a vendored HyperLogLog implementation.
- Record the software version or Git commit used for manuscript analyses.

Do not include here:

- HyperLogLog equations;
- minimizer definitions;
- the mixed-stride algorithm;
- benchmark-specific hardware or thread counts.

### 4.1.2 Input collections and pairwise execution

Include:

- The two positional inputs are text files containing paths to query and reference datasets.
- Hammock constructs one reusable sketch per listed file and emits a comma-separated table of pairwise comparisons between the two collections.
- Plain-text BED files are accepted by interval modes.
- FASTA inputs, including gzip-compressed FASTA, are accepted by sequence mode.
- BED inputs can be converted to FASTA by pairing each collection with a native reference through `--ref`, `--ref1`, or `--ref2`.
- Recognized compressed BED and BigBed inputs are rejected rather than silently interpreted as empty text.

The scientific point to retain is that file-level sketches are constructed once and reused across comparisons.

### 4.1.3 Implemented comparison modes

| CLI mode | Internal letter | Input | Representation |
|---|---:|---|---|
| `interval-string` | A | BED | Exact chromosome-start-end interval strings |
| `interval` / `interval-points` | B | BED | Base-level genomic positions covered by intervals |
| `interval-hybrid` | C | BED | Combined interval-string and position representation |
| `sequence` | D | FASTA or reference-backed BED | Sliding-window minimizers derived from sequence |

Manuscript emphasis:

- Mode B is the primary within-reference method.
- Mode D is the primary cross-reference method.
- Mixed-stride point subsampling applies to the point component of interval-points and interval-hybrid modes.
- Modes A and C should be acknowledged as implemented but not treated as headline contributions unless directly evaluated.

### 4.1.4 Reference-backed sequence conversion

Include:

- A BED collection can be paired with one shared reference or two collection-specific references.
- Sequence extraction is performed with `bedtools getfasta` against the reference associated with each BED collection.
- References may be supplied as local FASTA paths or resolved from a prepopulated cache.
- Reference downloads are performed separately through `hammock fetch-ref`; analysis runs do not download references automatically.
- Reference FASTAs are indexed with `samtools faidx` where required.
- Reference identifiers are retained in sequence-mode output.
- Cross-reference sequence similarity measures shared sequence content and should not be described as coordinate overlap, liftover, or proof of homologous loci.

The scientific details of sequence representation belong in Section 4.3.

### 4.1.5 Reproducibility and execution controls

Include:

- `--seed` controls xxh64 hashing for HLL ingestion.
- `--gate-seed` separately controls interval-point sampling and mixed-stride phase assignment.
- Mixed-stride output is deterministic for fixed input, `subB`, and gate seed.
- Thread count is user configurable.
- Output filenames encode the principal sketch and mode parameters.
- Exact benchmark thread counts, repetitions, hardware, and software environments belong in Section 4.6.

### 4.1.6 Output statistics

Include:

- `jaccard_similarity_ie`: set-Jaccard estimate obtained from HLL cardinality estimates using inclusion-exclusion. This is the principal interval-mode metric used for comparison with `bedtools jaccard`.
- `jaccard_similarity`: register-equality statistic, defined as the fraction of active HLL registers whose stored values are equal. It is retained for compatibility with the original Hammock implementation and is not ordinary set Jaccard.
- `containment_AB` and `containment_BA`: directional containment estimates.
- `cosketch_geom`, `cosketch_arith`, and `cosketch_max`: symmetric summaries of the directional containments.
- Sequence-mode rows record the reference associated with each collection when sequence is extracted from BED input.

Formal equations and censoring of negative inclusion-exclusion intersection estimates belong in Section 4.4.

## 4.2 Interval-mode sketch construction

### 4.2.1 Covered-position representation

- Define the element universe as base-level genomic positions covered by BED intervals.
- Specify how overlapping intervals within a file are handled.
- Define chromosome and coordinate encoding before hashing.
- Describe HLL precision and sketch size.

### 4.2.2 Mixed-stride point subsampling

- Define `subB` as the fraction of covered positions admitted to the sketch.
- Describe the mixed-stride rule, including stride length and chromosome-specific phase assignment.
- Explain why mixed stride avoids a per-position decision hash and reduces work approximately in proportion to the sampling rate.
- State determinism for fixed input, `subB`, and gate seed.
- Describe hash-threshold and single-hash samplers only as comparison methods used in evaluation.

## 4.3 Sequence-mode sketch construction

Proposed content:

- BED-to-FASTA extraction against each collection's native reference.
- FASTA record handling.
- Sliding-window minimizer selection.
- Definitions of `k` and `w`.
- Hashing and HLL ingestion of minimizers.
- One file-level sketch per interval-derived sequence collection.
- Interpretation as sequence-content similarity rather than coordinate overlap.

## 4.4 Similarity estimators and computational complexity

Proposed content:

- HLL cardinality estimation.
- Union by register-wise maximum.
- Inclusion-exclusion estimate of intersection and set Jaccard.
- Truncation of negative estimated intersections at zero.
- Register-equality statistic as a separate compatibility metric.
- Directional containment and cosketch summaries.
- Complexity of sketch construction and fixed-size sketch comparison.

## 4.5 Datasets

To be developed from the datasets actually used in Figures 3–7.

## 4.6 Performance benchmarking

To include hardware, software versions, compiler, thread counts, repetitions, timing procedure, and the exact mixed-stride settings used in Figure 3.

## 4.7 Interval-mode accuracy evaluation

Primary evaluation:

- Compare `jaccard_similarity_ie` with exact `bedtools jaccard`.
- Use identity-line agreement, Pearson correlation, rank correlation, and error summaries as appropriate.
- Register equality may be reported in supplementary analysis if useful, but should not compete with the principal set-Jaccard result in the main figure.

## 4.8 Sequence-mode biological evaluation

To include the cross-reference tissue experiment, Maurano tissue clustering, parameter sweeps, clustering procedure, and validation metrics.

# Figure 4 working decision

Recommended main-text design:

- Show inclusion–exclusion Jaccard against BEDTools Jaccard.
- Retain the identity line and make the caption explicit that Hammock values are `jaccard_similarity_ie` estimates derived from HLL cardinalities by inclusion-exclusion.
- Remove register equality from the main figure unless the paper needs a secondary panel specifically to explain backward compatibility.
- If register equality remains scientifically useful, move it to a supplementary figure or a short Methods note.

Rationale: Figure 4 should answer one clear question—whether Hammock accurately estimates the exact set-Jaccard quantity reported by BEDTools. Showing a second, intentionally different statistic complicates that claim without strengthening it.