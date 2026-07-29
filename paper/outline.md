## Thesis
Modern interval collections are limited by two boundaries: the computational cost of comparing every dataset and the coordinate dependence that prevents comparison across references. Hammock creates reusable sketches of either interval coordinates or interval-derived sequences, expanding comparison within and across those boundaries while preserving useful biological structure.

## Abstract

## I. Introduction

## II. Results
   ### 2.1 Hammock represents interval sets in complementary coordinate and sequence spaces

![Figure 2](figures/hammock_workflow_07292026.png)

**Figure 2. Overview of the hammock workflow.** (A) Public interval collections are large, comparisons scale quadratically with the number of files, and BED coordinates are tied to specific reference genomes. (B) In interval mode, each BED file is summarized as a sketch of covered genomic positions, enabling fast all-pairs similarity comparisons within a shared reference. (C) In sequence mode, interval sequences are extracted from each file’s native reference FASTA, summarized through minimizer-based sketches, and compared across references without requiring direct coordinate overlap.

   2.2 Interval sketching expands feasible all-pairs BED comparison
   2.3 Interval sketches preserve exact-overlap similarity structure
   2.4 Sequence sketches enable comparison across genome references
   2.5 Biological identity is preserved across references
   2.6 Sequence sketches recover tissue organization
   2.7 Parameterization separates numerical agreement from biological resolution

## III. Discussion
   3.1 Scaling comparison within references
   3.2 Comparing interval-derived sequence across references
   3.3 Coordinate similarity and sequence similarity are complementary
   3.4 Practical recommendations
   3.5 Limitations and future work

## IV. Methods
   4.1 Software implementation
   4.2 Interval-mode sketch construction
   4.3 Sequence-mode sketch construction
   4.4 Similarity estimators and computational complexity
   4.5 Datasets
   4.6 Performance benchmarking
   4.7 Interval-mode accuracy evaluation
   4.8 Sequence-mode biological evaluation

Data and code availability
Author contributions
Acknowledgments
References

Supplementary Methods
S1. Quantification of public interval collections
    S1.1 ChIP-Atlas experiment counts
    S1.2 ENCODE assay counts
    S1.3 Roadmap and BLUEPRINT file selection
    S1.4 Reference-assembly classification
    S1.5 Pairwise-comparison calculations
Supplementary Results
Supplementary Figures and Tables