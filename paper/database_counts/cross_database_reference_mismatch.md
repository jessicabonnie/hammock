# Cross-database reference mismatch case study

## Goal

Identify a biologically coherent collection of BED files that a researcher may want to compare, but cannot compare directly with coordinate-overlap tools because the files come from different reference assemblies.

## Primary case: H3K27ac in Roadmap Epigenomics and ENCODE

This case deliberately uses two public resources rather than two representations of the same experiment.

- **Roadmap Epigenomics** uniformly reprocessed human histone ChIP-seq data against **hg19** and distributes narrowPeak, broadPeak, and gappedPeak interval files.
- **ENCODE** currently distributes released H3K27ac peak BED files primarily on **GRCh38**.
- Both resources contain H3K27ac profiles across human tissues, primary cells, and cell lines.

The resulting files have the same broad biological purpose—identifying active regulatory regions—but do not share a coordinate system. A direct all-against-all BED overlap analysis across the two resources is therefore undefined unless one collection is lifted or remapped.

## Counting units

Report the following separately:

1. Roadmap H3K27ac peak BED files on hg19.
2. ENCODE H3K27ac peak BED files on GRCh38.
3. Distinct Roadmap epigenomes or samples represented.
4. Distinct ENCODE experiments represented.
5. Within-resource pair counts, which are directly comparable because each resource uses one assembly.
6. Cross-resource pairs, computed as `n_roadmap * n_encode`, which are biologically plausible comparisons but are not directly coordinate-comparable.

The cross-resource pair count is not a claim that every pair is equally biologically informative. It quantifies the number of candidate comparisons blocked by incompatible coordinate systems before any tissue or cell-type matching is applied.

## Peak representation

Use one peak representation consistently within each analysis. The first implementation uses narrowPeak files because both resources expose narrow peak calls. Broad and gapped H3K27ac representations should be counted separately rather than mixed into the same comparison set.

## Recommended refinement

After the repository-level count, create a stricter tissue- or cell-class-matched subset. For example:

- liver and hepatocyte H3K27ac;
- heart and cardiac-cell H3K27ac;
- brain-region H3K27ac;
- pluripotent stem-cell H3K27ac.

This refinement will show a situation where a user has a clear biological reason to compare files from both resources but conventional BED overlap still requires coordinate conversion.

## Candidate motivating statements

> Cross-database integration can fail even when every file describes the same assay target in the same species: Roadmap H3K27ac peaks are distributed on hg19, whereas current ENCODE peak files are distributed on GRCh38.

> Within each repository, H3K27ac BED files can be compared directly. Across repositories, the same analysis is blocked by incompatible reference coordinates unless one collection is lifted, remapped, or represented in a reference-independent form.

> The practical reference problem is not that a user would compare duplicate hg19 and hg38 versions of one experiment. It is that biologically related datasets from different resources may exist only on different assemblies.

> A reference-independent similarity representation would allow researchers to search across both collections without first forcing every BED file onto a single coordinate system.

## Sources

- Roadmap processed-data documentation states that Release 9 reads were mapped to hg19 and links the consolidated peak directories.
- ENCODE file records expose an explicit `assembly` field and provide released H3K27ac BED peak files on GRCh38.
- The NCBI Roadmap listing reports 164 H3K27ac records, but the final BED-file count must come from the specific Roadmap peak directory and selected peak representation rather than from the general experiment listing.
