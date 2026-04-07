# ChIP-seq reference comparison experiments

This directory holds ancillary data and workflows for ChIP-seq–related comparisons in Hammock. Layout:

```text
experiments/chip-seq-reference-comparison/
├── README.md                 # this file
├── data/
│   ├── histone-marks/        # ENCODE histone ChIP-seq cohort (narrowPeak + FASTQs + fasta_from_bed)
│   └── raw_files1.txt        # legacy raw URL list
├── manifest.txt
└── scripts/                  # helpers for histone-marks → nf-core/chipseq
    ├── build_histone_fastq_manifest.py
    ├── build_histone_narrowpeak_bed_manifest.py
    ├── build_nfcore_chipseq_samplesheet.py
    ├── chipseq_nfcore_samplesheet.py
    ├── histone_bedpaths_to_fasta_manifest.py  # manifest-driven bedtools getfasta → fasta_from_bed/
    ├── plot_hammock_histone_similarity.py    # hammock CSV + manifest → heatmaps / violins
    └── histone_encode_manifest_common.py   # shared cart parsing + batch columns
```

On some systems the same tree is mirrored under a shared filesystem, e.g. `/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/`. Commands below use a variable so you can point at either copy.

---

## Data preparation (histone marks, ENCODE)

All outputs are written under **`data/histone-marks/`** when you set `HISTONE_MARKS` as follows (from the **repository root** `hammock/`):

```bash
export HISTONE_MARKS="/home/jbonnie1/interval_sketch/hammock/experiments/chip-seq-reference-comparison/data/histone-marks"
# or, on a typical vast mirror:
# export HISTONE_MARKS="/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/data/histone-marks"
export HISTONE_MARKS_RESULTS="/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/results/histone-marks"
```

### 1. Cart export from ENCODE

1. Build an ENCODE cart (e.g. mouse/human histone ChIP with FASTQs and narrowPeak beds).
2. Download the **Experiment** cart report as TSV (tabular). Save it as:

   ```text
   ${HISTONE_MARKS}/experiment_report_2026_3_23_17h_57m.tsv
   ```

   (Use your actual export date in the filename if you regenerate it; the scripts only need a consistent path you pass on the command line.)

3. From the same project, export **processed** file URLs (e.g. `histone_bed.txt` for beds + bigWigs etc.) and **raw** FASTQ URLs (`histone_raw.txt`). Place them in `${HISTONE_MARKS}/`.

### 2. Download IP FASTQs

```bash
mkdir -p "${HISTONE_MARKS}/fastq"
# Skip the first line (metadata URL); download each remaining URL into fastq/
tail -n +2 "${HISTONE_MARKS}/histone_raw.txt" | while read -r url; do
  case "$url" in http*) curl -fL --retry 3 -o "${HISTONE_MARKS}/fastq/$(basename "$url")" "$url" ;; esac
done
```

### 3. Download narrowPeak BED files (optional for sketching / benchmarks)

Only **`.bed.gz`** narrowPeaks are needed if you are not using bigWig/BAM yet:

```bash
mkdir -p "${HISTONE_MARKS}/bed"
grep '\.bed\.gz$' "${HISTONE_MARKS}/histone_bed.txt" | while read -r url; do
  curl -fL --retry 3 -o "${HISTONE_MARKS}/bed/$(basename "$url")" "$url"
done
```

### 4. narrowPeak BED manifest (TSV linking `bed/*.bed.gz` to cart metadata)

From the repository root. Writes **`${HISTONE_MARKS}/histone_narrowpeak_bed_manifest.tsv`** (one row per local `ENCFFxxxxxx.bed.gz` that appears in the cart **Files** column):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_histone_narrowpeak_bed_manifest.py \
  --experiment-report "${HISTONE_MARKS}/experiment_report_2026_3_23_17h_57m.tsv" \
  --bed-dir "${HISTONE_MARKS}/bed" \
  --output "${HISTONE_MARKS}/histone_narrowpeak_bed_manifest.tsv"
```

**With ENCODE per-file assembly** (slow: one HTTPS GET per local BED). Use this when the cart mixes hg19/GRCh38 or mm9/mm10/mm10-minimal—then run **`histone_bedpaths_to_fasta_manifest.py --manifest`** (§5):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_histone_narrowpeak_bed_manifest.py \
  --experiment-report "${HISTONE_MARKS}/experiment_report_2026_3_23_17h_57m.tsv" \
  --bed-dir "${HISTONE_MARKS}/bed" \
  --output "${HISTONE_MARKS}/histone_narrowpeak_bed_manifest.tsv" \
  --fetch-assembly
```

With **`--fetch-assembly`**, the **`assembly`** column is filled from ENCODE (`mm10`, `mm10-minimal`, `mm9`, `GRCh38`, `hg19`, …) and **`reference_build`** is set to the same token. Without it, **`assembly`** stays empty and **`reference_build`** defaults to **mm10** / **GRCh38** from cart **`organism`** only (wrong for mixed-assembly carts).

Columns include `bed_path`, `file_accession`, **`assembly`**, **`reference_build`**, `experiment_accession`, batch/provenance fields (`lab`, `project`, `library_construction_platform`, `library_construction_method` from the cart), biosample fields, `tissue_type`, `histone_target`, `organism`, replicate fields, and `linked_antibody`.

The same run also writes **one BED path per line per reference build** under **`${HISTONE_MARKS}/path_lists/`** (for tools that want plain list files), e.g.:

- **`path_lists/histone_narrowpeak_bed_paths_mm10.txt`**, **`path_lists/histone_narrowpeak_bed_paths_GRCh38.txt`**, **`path_lists/histone_narrowpeak_bed_paths_hg19.txt`**, **`path_lists/histone_narrowpeak_bed_paths_mm9.txt`**, **`path_lists/histone_narrowpeak_bed_paths_mm10_minimal.txt`**

Only builds that appear in the manifest (after resolving `assembly` / `reference_build` / `organism`) get a file. Rows with no resolvable build are omitted (a warning reports the count).

**Re-run the manifest builder** when you add or remove local BEDs or refresh the cart TSV. Add **`--fetch-assembly`** when you need ENCODE’s true per-file assembly (mixed hg19/GRCh38, mm9/mm10/mm10-minimal); without it, `assembly` stays blank and lists are grouped only by the organism fallback (**mm10** / **GRCh38**).

### 4b. broadPeak BED path lists (by build)

If your local `bed/` directory also contains broadPeak `.bed.gz` files, build per-build broadPeak lists from `bed/manifest.tsv`:

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_histone_broadpeak_bed_paths.py \
  --bed-manifest "${HISTONE_MARKS}/bed/manifest.tsv" \
  --bed-dir "${HISTONE_MARKS}/bed" \
  --path-lists-dir "${HISTONE_MARKS}/path_lists"
```

### 5. FASTA from narrowPeak BEDs (`fasta_from_bed/`)

Create **`${HISTONE_MARKS}/fasta_from_bed/`** (sibling of `bed/` and `fastq/`) with one **`ENCFFxxxxxx.fa`** per narrowPeak **`ENCFFxxxxxx.bed.gz`**, using **`bedtools getfasta`** and a reference that matches each file’s ENCODE **`assembly`**.

**Mixed assemblies.** A single histone cart often includes **both hg19 and GRCh38** human peaks and **mm9, mm10, and mm10-minimal** mouse peaks. Using only GRCh38 + mm10 for every file mis-sequences a large fraction of peaks and can produce `bedtools` warnings (intervals past the end of a chromosome in your `-fi` FASTA). Prefer **`--manifest`** with a manifest built using **`--fetch-assembly`** (§4), and set FASTAs for every assembly that appears.

| ENCODE `assembly` | Default FASTA under `/data/blangme2/jessica/mus_homo/references/` | Override |
|-------------------|------------------------------------------------------------------------|----------|
| **GRCh38** | `grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` (symlink) | `$GRCH38_FASTA` / `--ref-grch38` |
| **hg19** | `hg19/hg19.fa` | `$HG19_FASTA` / `--ref-hg19` |
| **mm10** | `mm10/mm10_no_alt_analysis_set_ENCODE.fasta` (via `references/mm10` → `…/jessica/mm10`) | `$MM10_FASTA` / `--ref-mm10` |
| **mm10-minimal** | `mm10_minimal/mm10_minimal_ENCODE.fasta` | `$MM10_MINIMAL_FASTA` / `--ref-mm10-minimal` |
| **mm9** | `mm9/mm9.fa` | `$MM9_FASTA` / `--ref-mm9` |

Older mus-homo docs also mention [`experiments/mus-homo/dnase-seq/dnase-seq.md`](mus-homo/dnase-seq/dnase-seq.md) paths; the table above is what **`histone_bedpaths_to_fasta_manifest.py`** uses by default.

Requires **`bedtools`** on `PATH`. Gzipped BED is passed through to `bedtools` as in the mus-homo recipes.

From the **repository root** (`hammock/`), after **`HISTONE_MARKS`** is set:

**Recommended (assembly-aware):** rebuild the manifest with **`--fetch-assembly`** if needed, then (defaults use the **`mus_homo/references/`** tree above):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/histone_bedpaths_to_fasta_manifest.py \
  --manifest "${HISTONE_MARKS}/histone_narrowpeak_bed_manifest.tsv" \
  --histone-marks-dir "${HISTONE_MARKS}"
```

Add **`--ref-hg19`**, **`--ref-mm9`**, etc., only if those files live somewhere else on your system.

If the **`assembly`** column is missing or entirely blank, **`histone_bedpaths_to_fasta_manifest.py`** still resolves references from **`reference_build`** (and finally from **`organism`**). A warning is printed when **`assembly`** was never filled from ENCODE (**`--fetch-assembly`**).

The script writes FASTAs under **`${HISTONE_MARKS}/fasta_from_bed/narrowpeak/`** and **`${HISTONE_MARKS}/fasta_from_bed/broadpeak/`** (separate directories to avoid collisions), and refreshes **`path_lists/histone_*_fasta_paths_<build>.txt`**. In manifest mode it writes **both** species-based narrowPeak lists (**`path_lists/histone_narrowpeak_fasta_paths_Mus_musculus.txt`**, **`_Homo_sapiens.txt`**) and per-build narrowPeak lists, and it also processes broadPeak build lists from **`path_lists/histone_broadpeak_bed_paths_<build>.txt`** when present. It also writes **`tissue_type_breakout_by_build_and_organism.png`**. Use **`--skip-existing`**, **`--dry-run`**, and **`-j N`** as needed (parallel `bedtools`; **`-j 1`** is sequential). Disable PNG generation with **`--no-tissue-plot`**.

### 6. FASTQ manifest (TSV linking paths to cart metadata)

From the repository root (so `scripts/` resolves correctly):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_histone_fastq_manifest.py \
  --experiment-report "${HISTONE_MARKS}/experiment_report_2026_3_23_17h_57m.tsv" \
  --fastq-dir "${HISTONE_MARKS}/fastq" \
  --output "${HISTONE_MARKS}/histone_fastq_manifest.tsv"
```

This writes **`histone_fastq_manifest.tsv`** next to the other histone-marks files (same column set as the bed manifest: paths, accessions, **`lab` / `project` / library construction** fields for batch awareness, biosample and replicate fields, `linked_antibody`).

### 7. nf-core/chipseq samplesheets + control downloads (mouse vs human)

Requires outbound HTTPS to `encodeproject.org` (throttled by `--sleep`).

Samplesheets are **split by species**: **Mus musculus** and **Homo sapiens** get separate CSVs, control URL tables, and download scripts so you can run [nf-core/chipseq](https://nf-co.re/chipseq/usage#samplesheet-input) with one `--genome` per workflow (e.g. mm10 vs GRCh38).

Use the FASTQ manifest (§6) so the `organism` column is available; if you omit `--manifest`, organism is inferred from each IP file’s experiment `biosample_summary` on ENCODE (extra API calls).

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_nfcore_chipseq_samplesheet.py \
  --fastq-dir "${HISTONE_MARKS}/fastq" \
  --output-dir "${HISTONE_MARKS}" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv"
```

This creates, under **`${HISTONE_MARKS}/`**:

| Output | Purpose |
|--------|---------|
| `encode_file_metadata_cache.jsonl` | Cached ENCODE `/files/ENCFF…/` JSON |
| `nfcore_chipseq_samplesheet_Mus_musculus.csv` | Mouse-only `--input` |
| `nfcore_chipseq_samplesheet_Homo_sapiens.csv` | Human-only `--input` (omitted if no human FASTQs) |
| `nfcore_control_fastqs_encode_Mus_musculus.tsv` | Mouse input FASTQ URLs + paths |
| `nfcore_control_fastqs_encode_Homo_sapiens.tsv` | Human input FASTQ URLs + paths |
| `download_control_fastqs_Mus_musculus.sh` | `curl` mouse controls into `fastq_controls/` |
| `download_control_fastqs_Homo_sapiens.sh` | `curl` human controls |

Each CSV includes an extra column **`encode_organism`** (ignored by nf-core; useful for checks).

Run-specific note (2026-04-07): after generating these samplesheets, the mouse file was manually edited to drop failing replicate 2 for `ENCSR548BKP_H3K4me1` (rows with IP FASTQs `ENCFF089POY` and `ENCFF980EMT`). This replicate failed in `nf-core/chipseq` during MACS3 model-building (`<100` paired peaks); replicate 1 for the same experiment was retained.

Run-specific note (2026-04-07, later rerun): the mouse file was edited again to drop failing replicate 1 for `ENCSR282GAG_H3K27me3` (rows with IP FASTQs `ENCFF213QFH` and `ENCFF683GIH`) after another MACS3 model-building failure (`94` paired peaks, still below the `100` minimum). Replicate 2 for this experiment was retained.

Fetch controls (default destination: **`${HISTONE_MARKS}/fastq_controls/`**):

```bash
bash "${HISTONE_MARKS}/download_control_fastqs_Mus_musculus.sh"
bash "${HISTONE_MARKS}/download_control_fastqs_Homo_sapiens.sh"   # if present
```

<!-- Then run Nextflow **per species**, for example:

```bash
nextflow run chipseq -r 2.1.0 \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Mus_musculus.csv" \
  --outdir "${HISTONE_MARKS_RESULTS}/nf-core_mm10" \
  -profile singularity \
  --genome mm10 \
  --macs_gsize 1870000000

nextflow run chipseq -r 2.1.0 \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --outdir "${HISTONE_MARKS_RESULTS}/nf-core-GRCh38" \
  -profile singularity \
  --genome GRCh38 \
  --macs_gsize 2700000000

  nextflow run nf-core/chipseq -r 2.1.0 \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --outdir "${HISTONE_MARKS_RESULTS}/nf-core-GRCh37" \
  -profile singularity \
  --genome GRCh37 \
  --macs_gsize 2700000000
``` -->

Optional **single combined** sheet (all species in one file, previous behavior):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_nfcore_chipseq_samplesheet.py \
  --fastq-dir "${HISTONE_MARKS}/fastq" \
  --output-dir "${HISTONE_MARKS}" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv" \
  --write-combined
```

Re-runs with a warm cache:

```bash
python3 experiments/chip-seq-reference-comparison/scripts/build_nfcore_chipseq_samplesheet.py \
  --fastq-dir "${HISTONE_MARKS}/fastq" \
  --output-dir "${HISTONE_MARKS}" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv" \
  --skip-fetch
```

(`--skip-fetch` avoids re-downloading IP file metadata from ENCODE if the JSONL cache is already complete; organism resolution and control pairing may still contact ENCODE when needed.)

### Summary of manifest paths under `${HISTONE_MARKS}/`

| File | Produced by |
|------|-------------|
| `histone_narrowpeak_bed_manifest.tsv` | `build_histone_narrowpeak_bed_manifest.py` (§4) |
| `path_lists/histone_narrowpeak_bed_paths_<build>.txt` | §4 (per reference build, e.g. `mm10`, `GRCh38`) |
| `path_lists/histone_broadpeak_bed_paths_<build>.txt` | §4b (per reference build) |
| `fasta_from_bed/narrowpeak/*.fa` + `fasta_from_bed/broadpeak/*.fa` | `histone_bedpaths_to_fasta_manifest.py` (§5) |
| `path_lists/histone_*_fasta_paths_<build>.txt` | §5 (`--manifest`; same slug as BED list) |
| `path_lists/histone_narrowpeak_fasta_paths_Mus_musculus.txt` / `_Homo_sapiens.txt` | §5 (**`--manifest`** only; species) |
| `histone_fastq_manifest.tsv` | `build_histone_fastq_manifest.py` (§6) |

---

## Scripts reference

| Script | Role |
|--------|------|
| `scripts/histone_encode_manifest_common.py` | Shared: cart parsing, batch fields, `normalize_organism_label`, species BED path list writer |
| `scripts/build_histone_narrowpeak_bed_manifest.py` | Cart TSV + `bed/*.bed.gz` → `histone_narrowpeak_bed_manifest.tsv` |
| `scripts/build_histone_fastq_manifest.py` | Cart TSV + `fastq/` → `histone_fastq_manifest.tsv` |
| `scripts/histone_bedpaths_to_fasta_manifest.py` | Manifest + assembly-aware refs → `fasta_from_bed/{narrowpeak,broadpeak}/*.fa` + FASTA path lists |
| `scripts/build_nfcore_chipseq_samplesheet.py` | Local IP FASTQs + manifest → per-species nf-core CSVs + control TSVs + downloaders |
| `scripts/chipseq_nfcore_samplesheet.py` | Library: ENCODE API, samplesheet build, split by `Mus musculus` / `Homo sapiens` |
| `scripts/plot_hammock_histone_similarity.py` | Hammock pairwise CSV + `histone_fastq_manifest.tsv` (+ optional nf-core samplesheet) → annotated TSV, clustermaps, tissue blocks, cross-reference violins |

---

## Reference-genome comparison workflow 

The following describes an older **same-TF, two-reference** alignment comparison (not the histone-marks ENCODE cohort). Placeholders remain for cell line, TF, and references.

### Objective

Compare ChIP-seq data for the same cell line and transcription factor aligned to different reference genomes using Hammock Mode D.

### Experimental design

- **Transcription factor**: [TO BE FILLED]
- **Reference genome 1**: [TO BE FILLED]
- **Reference genome 2**: [TO BE FILLED]
- **Number of replicates**: [TO BE FILLED]

### Data treatment

#### 1. Data preparation

##### Run nf-core/chipseq for both references

```bash
nextflow run chipseq/ \
  -profile singularity \
  --outdir ${HISTONE_MARKS_RESULTS}/nf-core-GRCh37/ \
  --input ${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv \
  --genome GRCh37 \
  --macs_gsize 2700000000

nextflow run chipseq/ \
  -profile singularity \
  --outdir ${HISTONE_MARKS_RESULTS}/nf-core-GRCh38/ \
  --input ${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv \
  --genome GRCh38 \
  --macs_gsize 2700000000

nextflow run chipseq \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Mus_musculus.csv" \
  --outdir "${HISTONE_MARKS_RESULTS}/nf-core-mm10" \
  -profile singularity \
  --genome mm10 \
  --macs_gsize 1870000000
```

#### 2. Organize BED files

```bash
# cd chip-seq-reference-comparison/results/histone-marks

# Adjust these to your actual nf-core output directories
REF_G37_OUT=${HISTONE_MARKS_RESULTS}/nf-core-GRCh37
REF_G38_OUT=${HISTONE_MARKS_RESULTS}/nf-core-GRCh38/
REF_MM10_OUT=${HISTONE_MARKS_RESULTS}/nf-core-mm10

# Recursive search (bash ** only works with shopt -s globstar; nf-core nests peaks several dirs deep).
# xargs -r avoids "realpath: missing operand" when find returns nothing (pipeline not finished or wrong --outdir).
find "$REF_G37_OUT" -type f -name '*_peaks.broadPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_GRCh37_bed_paths_BROAD.txt
find "$REF_G37_OUT" -type f -name '*_peaks.narrowPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_GRCh37_bed_paths_NARROW.txt

find "$REF_G38_OUT" -type f -name '*_peaks.broadPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_GRCh38_bed_paths_BROAD.txt
find "$REF_G38_OUT" -type f -name '*_peaks.narrowPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_GRCh38_bed_paths_NARROW.txt

find "$REF_MM10_OUT" -type f -name '*_peaks.broadPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_mm10_bed_paths_BROAD.txt
find "$REF_MM10_OUT" -type f -name '*_peaks.narrowPeak' -print0 2>/dev/null | xargs -0r realpath > ${HISTONE_MARKS_RESULTS}/ref_mm10_bed_paths_NARROW.txt

```

#### 3. Create FASTA files from peaks

```bash
mkdir -p ${REF_G38_OUT}/fastas_from_broad ${REF_G38_OUT}/fastas_from_narrow
mkdir -p ${REF_G37_OUT}/fastas_from_broad ${REF_G37_OUT}/fastas_from_narrow
mkdir -p ${REF_MM10_OUT}/fastas_from_broad ${REF_MM10_OUT}/fastas_from_narrow

python3 experiments/chip-seq-reference-comparison/scripts/beds_to_fastas.py \
  --bed-list ${HISTONE_MARKS_RESULTS}/ref_GRCh37_bed_paths_BROAD.txt \
  --output-dir ${REF_G37_OUT}/fastas_from_broad \
  --reference GRCh37 \
  -j 8 > ${HISTONE_MARKS_RESULTS}/ref_GRCh37_fasta_paths_BROAD.txt
  # use dry run if only the fasta paths are needed as output

  python3 experiments/chip-seq-reference-comparison/scripts/beds_to_fastas.py \
  --bed-list ${HISTONE_MARKS_RESULTS}/ref_GRCh38_bed_paths_BROAD.txt \
  --output-dir ${REF_G38_OUT}/fastas_from_broad \
  --reference GRCh38 \
  -j 8 > ${HISTONE_MARKS_RESULTS}/ref_GRCh38_fasta_paths_BROAD.txt

  python3 experiments/chip-seq-reference-comparison/scripts/beds_to_fastas.py \
  --bed-list ${HISTONE_MARKS_RESULTS}/ref_GRCh37_bed_paths_NARROW.txt \
  --output-dir ${REF_G37_OUT}/fastas_from_narrow \
  --reference GRCh37 \
  -j 8 > ${HISTONE_MARKS_RESULTS}/ref_GRCh37_fasta_paths_NARROW.txt

  python3 experiments/chip-seq-reference-comparison/scripts/beds_to_fastas.py \
  --bed-list ${HISTONE_MARKS_RESULTS}/ref_GRCh38_bed_paths_NARROW.txt \
  --output-dir ${REF_G38_OUT}/fastas_from_narrow \
  --reference GRCh38 \
  -j 8 > ${HISTONE_MARKS_RESULTS}/ref_GRCh38_fasta_paths_NARROW.txt

  python3 experiments/chip-seq-reference-comparison/scripts/beds_to_fastas.py \
  --bed-list ${HISTONE_MARKS_RESULTS}/ref_mm10_bed_paths_BROAD.txt \
  --output-dir ${REF_MM10_OUT}/fastas_from_broad \
  --reference mm10 \
  -j 8 > ${HISTONE_MARKS_RESULTS}/ref_mm10_fasta_paths_BROAD.txt

cat ${HISTONE_MARKS_RESULTS}/ref_GRCh38_fasta_paths_BROAD.txt ${HISTONE_MARKS_RESULTS}/ref_GRCh37_fasta_paths_BROAD.txt > ${HISTONE_MARKS_RESULTS}/homo_BROAD_fastas.txt

cat ${HISTONE_MARKS_RESULTS}/ref_GRCh38_fasta_paths_NARROW.txt ${HISTONE_MARKS_RESULTS}/ref_GRCh37_fasta_paths_NARROW.txt > ${HISTONE_MARKS_RESULTS}/homo_NARROW_fastas.txt

cat ${HISTONE_MARKS_RESULTS}/homo_BROAD_fastas.txt ${HISTONE_MARKS_RESULTS}/ref_mm10_fasta_paths_BROAD.txt > ${HISTONE_MARKS_RESULTS}/mus_homo_BROAD_fastas.txt

```

### Mode D analysis

#### Humans Across References -- BROAD/NARROW Peak
```bash
mkdir ${HISTONE_MARKS_RESULTS}/hammock/

hammock ${HISTONE_MARKS_RESULTS}/homo_BROAD_fastas.txt \
  ${HISTONE_MARKS_RESULTS}/homo_BROAD_fastas.txt \
    --outprefix ${HISTONE_MARKS_RESULTS}/hammock/homo_BROAD \
    --minimizer \
    -k 5 -w 40 --precision 24 \
    --full-paths

hammock ${HISTONE_MARKS_RESULTS}/homo_NARROW_fastas.txt \
  ${HISTONE_MARKS_RESULTS}/homo_NARROW_fastas.txt \
  --outprefix ${HISTONE_MARKS_RESULTS}/hammock/homo_NARROW \
  --minimizer \
  -k 5 -w 40 --precision 24 \
  --full-paths
```

The emitted CSV name includes sketch parameters (e.g. ``..._mnmzr_p24_jaccD_k5_w40.csv`` for minimizer runs with ``--precision 24``, ``-k 5``, ``-w 40``). Use **`--full-paths`** so paths join to the manifest.

#### Plotting (manifest-joined similarity)

**Biological replicates.** By default the plotting script uses **all** replicates (`--biological-replicate all`). Pairwise similarities are **averaged** when multiple file pairs correspond to the same biosample unit (same experiment, histone, tissue, and life stage, within a reference or across references for violins). Pass e.g. `--biological-replicate 1` to restrict to one ENCODE biological replicate after joining paths to metadata.

**Path filter.** Use **`--path-filter fastas_from_broad`** for broadPeak FASTAs (`.../fastas_from_broad/...`) or **`--path-filter fastas_from_narrow`** for narrowPeak FASTAs (`.../fastas_from_narrow/...`). The substring must appear in both `file1` and `file2` paths.

**Similarity column.** The plotting script always reads **`jaccard_similarity_with_ends`** (minimizer Jaccard including end k-mers). The CSV must contain that column, or pass **`--similarity-col`** to use a different column.

**2D embeddings (MDS + PCA).** Scatter plots share one scheme: **face color = tissue**, **marker shape = histone mark**, **edge color = reference build**, and **semi-transparent points** so overlaps are visible. **MDS:** **`mds_embedding_by_lifestage.pdf`** (one page per life stage; classical MDS on distance = 1 − similarity for that stage’s track matrix), **`mds_coordinates_by_lifestage.tsv`**, and **`similarity_matrix_all_tracks.tsv`** (full aggregated similarity matrix for reference). **PCA:** **`pca_by_lifestage.pdf`** (one page per life stage; PCA on centered **rows** of the life-stage similarity matrix), plus **`pca_coordinates_by_lifestage.tsv`**. Same similarity column as the rest of the pipeline.

**Outlier experiments.** In this histone reference dataset, two accessions were extreme on adult PCA after reference-aware sketching: **ENCSR330MAM** (H3K27me3, stomach, adult) and **ENCSR158WBG** (H3K4me1, stomach, adult; especially vs other experiments). By default the plotting script drops any pairwise row involving those **`experiment_accession`** values via **`--exclude-experiment-accession`** (default `ENCSR330MAM,ENCSR158WBG`). Pass **`--exclude-experiment-accession -`** to keep every experiment. The filter applies to all figures (not only PCA).

**Python environment.** The plot script needs `pandas`, `matplotlib`, `scipy`, and `seaborn`. If the system `python3` is too old or lacks those packages, use a conda env (example below) or `pip install pandas matplotlib seaborn scipy`.

**Commands — broadPeak vs narrowPeak.** From the repository root (`hammock/`), set paths and run **one** of the two blocks below. They differ only in the input CSV (`homo_BROAD_...` vs `homo_NARROW_...`), **`--path-filter`** (`fastas_from_broad` vs `fastas_from_narrow`), and **`--outdir`** (`figs-broad` vs `figs-narrow`). Adjust the CSV basename if your minimizer run used different `p`, `k`, or `w`.

```bash
# From repository root hammock/
export HISTONE_MARKS="/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/data/histone-marks"
export HISTONE_MARKS_RESULTS="/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/results/histone-marks"

mkdir -p "${HISTONE_MARKS_RESULTS}/hammock/figs-broad"

# broadPeak FASTAs (path segment .../fastas_from_broad/... in hammock CSV)
python3 experiments/chip-seq-reference-comparison/scripts/plot_hammock_histone_similarity.py \
  "${HISTONE_MARKS_RESULTS}/hammock/homo_BROAD_mnmzr_p24_jaccD_k5_w40.csv" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv" \
  --samplesheet "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --organism "Homo sapiens" \
  --path-filter fastas_from_broad \
  --outdir "${HISTONE_MARKS_RESULTS}/hammock/figs-broad"
```

```bash
mkdir -p "${HISTONE_MARKS_RESULTS}/hammock/figs-narrow"

# narrowPeak FASTAs (path segment .../fastas_from_narrow/...)
python3 experiments/chip-seq-reference-comparison/scripts/plot_hammock_histone_similarity.py \
  "${HISTONE_MARKS_RESULTS}/hammock/homo_NARROW_mnmzr_p24_jaccD_k5_w40.csv" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv" \
  --samplesheet "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --organism "Homo sapiens" \
  --path-filter fastas_from_narrow \
  --outdir "${HISTONE_MARKS_RESULTS}/hammock/figs-narrow"
```

The default is `--biological-replicate all` (mean across replicates as described above). Other outputs under `--outdir` include `annotated_pairs.tsv`, `cross_reference_violin.png` / `cross_reference_violin_faceted.png`, `clustermap_<ref>.png`, `clustermap_lifestage_<slug>.png` (life stage: tick labels show reference, histone, experiment; strips = tissue), `dendrogram_*.png`, `tissue_block_heatmap_*.png`, `similarity_matrix_*.tsv`, `tissue_mean_block_*.tsv`, and the MDS/PCA embedding files named above.

On a cluster where `python3` is 3.6 without pandas, use a newer interpreter, for example:

```bash
/data/apps/extern/anaconda/envs/mamba/0.23.0/bin/python3.10 \
  experiments/chip-seq-reference-comparison/scripts/plot_hammock_histone_similarity.py \
  "${HISTONE_MARKS_RESULTS}/hammock/homo_BROAD_mnmzr_p24_jaccD_k5_w40.csv" \
  --manifest "${HISTONE_MARKS}/histone_fastq_manifest.tsv" \
  --samplesheet "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --organism "Homo sapiens" \
  --path-filter fastas_from_broad \
  --outdir "${HISTONE_MARKS_RESULTS}/hammock/figs-broad"
```

