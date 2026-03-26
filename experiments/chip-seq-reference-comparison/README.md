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
    ├── histone_bedpaths_to_fasta.py       # bedtools getfasta → fasta_from_bed/
    └── histone_encode_manifest_common.py   # shared cart parsing + batch columns
```

On some systems the same tree is mirrored under a shared filesystem, e.g. `/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/`. Commands below use a variable so you can point at either copy.

---

## Data preparation (histone marks, ENCODE)

All outputs are written under **`data/histone-marks/`** when you set `HISTONE_MARKS` as follows (from the **repository root** `hammock/`):

```bash
export HISTONE_MARKS="$(pwd)/experiments/chip-seq-reference-comparison/data/histone-marks"
# or, on a typical vast mirror:
# export HISTONE_MARKS="/vast/blangme2/jbonnie/hammock/chip-seq-reference-comparison/data/histone-marks"
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

**With ENCODE per-file assembly** (slow: one HTTPS GET per local BED). Use this when the cart mixes hg19/GRCh38 or mm9/mm10/mm10-minimal—then run **`histone_bedpaths_to_fasta.py --manifest`** (§5):

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

Older mus-homo docs also mention [`experiments/mus-homo/dnase-seq/dnase-seq.md`](mus-homo/dnase-seq/dnase-seq.md) paths; the table above is what **`histone_bedpaths_to_fasta.py`** uses by default.

Requires **`bedtools`** on `PATH`. Gzipped BED is passed through to `bedtools` as in the mus-homo recipes.

From the **repository root** (`hammock/`), after **`HISTONE_MARKS`** is set:

**Recommended (assembly-aware):** rebuild the manifest with **`--fetch-assembly`** if needed, then (defaults use the **`mus_homo/references/`** tree above):

```bash
python3 experiments/chip-seq-reference-comparison/scripts/histone_bedpaths_to_fasta.py \
  --manifest "${HISTONE_MARKS}/histone_narrowpeak_bed_manifest.tsv" \
  --histone-marks-dir "${HISTONE_MARKS}"
```

Add **`--ref-hg19`**, **`--ref-mm9`**, etc., only if those files live somewhere else on your system.

If the **`assembly`** column is missing or entirely blank, **`histone_bedpaths_to_fasta.py`** still resolves references from **`reference_build`** (and finally from **`organism`**). A warning is printed when **`assembly`** was never filled from ENCODE (**`--fetch-assembly`**).

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

Fetch controls (default destination: **`${HISTONE_MARKS}/fastq_controls/`**):

```bash
bash "${HISTONE_MARKS}/download_control_fastqs_Mus_musculus.sh"
bash "${HISTONE_MARKS}/download_control_fastqs_Homo_sapiens.sh"   # if present
```

Then run Nextflow **per species**, for example:

```bash
nextflow run nf-core/chipseq -r 2.1.0 \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Mus_musculus.csv" \
  --outdir "${HISTONE_MARKS}/nfcore_chipseq_results_mouse" \
  -profile singularity \
  --genome mm10

nextflow run nf-core/chipseq -r 2.1.0 \
  --input "${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv" \
  --outdir "${HISTONE_MARKS}/nfcore_chipseq_results_human" \
  -profile singularity \
  --genome GRCh38
```

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
| `fasta_from_bed/*.fa` | `histone_bedpaths_to_fasta.py` (§5) |
| `path_lists/histone_*_fasta_paths_<build>.txt` | §5 (legacy list mode or **`--manifest`**; same slug as BED list) |
| `path_lists/histone_narrowpeak_fasta_paths_Mus_musculus.txt` / `_Homo_sapiens.txt` | §5 (**`--manifest`** only; species) |
| `histone_fastq_manifest.tsv` | `build_histone_fastq_manifest.py` (§6) |

---

## Scripts reference

| Script | Role |
|--------|------|
| `scripts/histone_encode_manifest_common.py` | Shared: cart parsing, batch fields, `normalize_organism_label`, species BED path list writer |
| `scripts/build_histone_narrowpeak_bed_manifest.py` | Cart TSV + `bed/*.bed.gz` → `histone_narrowpeak_bed_manifest.tsv` |
| `scripts/build_histone_fastq_manifest.py` | Cart TSV + `fastq/` → `histone_fastq_manifest.tsv` |
| `scripts/histone_bedpaths_to_fasta.py` | BED path lists + mm10/GRCh38 refs → `fasta_from_bed/*.fa` + FASTA path lists |
| `scripts/build_nfcore_chipseq_samplesheet.py` | Local IP FASTQs + manifest → per-species nf-core CSVs + control TSVs + downloaders |
| `scripts/chipseq_nfcore_samplesheet.py` | Library: ENCODE API, samplesheet build, split by `Mus musculus` / `Homo sapiens` |

---

## Reference-genome comparison workflow (separate design)

The following describes an older **same-TF, two-reference** alignment comparison (not the histone-marks ENCODE cohort). Placeholders remain for cell line, TF, and references.

### Objective

Compare ChIP-seq data for the same cell line and transcription factor aligned to different reference genomes using Hammock Mode D.

### Experimental design

- **Cell line**: [TO BE FILLED]
- **Transcription factor**: [TO BE FILLED]
- **Reference genome 1**: [TO BE FILLED]
- **Reference genome 2**: [TO BE FILLED]
- **Number of replicates**: [TO BE FILLED]

### Data treatment

#### 1. Data preparation

```bash
cd data
xargs -n 1 curl -O -L < raw_files1.txt
```

##### Download/locate reference genomes

```bash
# Reference 1 (e.g., hg19)
REF1=/path/to/reference1.fa

# Reference 2 (e.g., hg38)
REF2=/path/to/reference2.fa
```

##### Run nf-core/chipseq for both references

```bash
nextflow run chipseq/ -profile singularity --outdir chip-seq-reference-comparison/results/histone-marks/nf-core-GRCh37/ --input chip-seq-reference-comparison/data/histone-marks/nfcore_chipseq_samplesheet_Homo_sapiens.csv --genome GRCh37 --macs_gsize 2700000000
nextflow run chipseq/ -profile singularity --outdir chip-seq-reference-comparison/results/histone-marks/nf-core/ --input chip-seq-reference-comparison/data/histone-marks/nfcore_chipseq_samplesheet_Homo_sapiens.csv --genome GRCh38 --macs_gsize 2700000000
```

#### 2. Organize BED files

```bash
cd chip-seq-reference-comparison/results/histone-marks

# Adjust these to your actual nf-core output directories
REF_G37_OUT=nf-core-GRCh37
REF_G38_OUT=nf-core-GRCh38

# Typical nf-core/chipseq peak directories; update glob if your profile/version differs
ls "$REF_G37_OUT"/**/*_peaks.narrowPeak 2>/dev/null | xargs realpath > ref_GRCh37_bed_paths.txt
ls "$REF_G38_OUT"/**/*_peaks.narrowPeak 2>/dev/null | xargs realpath > ref_GRCh38_bed_paths.txt

while read -r filepath; do
    basename "$filepath" | sed 's/_peaks.narrowPeak$//'
done < ref_GRCh37_bed_paths.txt > sample_names.txt
```

#### 3. Create FASTA files from peaks

```bash
mkdir -p fastas_GRCh37 fastas_GRCh38

REF_G37=/path/to/GRCh37.fa
REF_G38=/path/to/GRCh38.fa

while read -r bedfile; do
    outfile=$(basename "$bedfile" _peaks.narrowPeak)
    bedtools getfasta -fi "$REF_G37" -bed "$bedfile" -fo "fastas_GRCh37/${outfile}.fa"
    echo "$(realpath "fastas_GRCh37/${outfile}.fa")"
done < ref_GRCh37_bed_paths.txt > ref_GRCh37_fastas.txt

while read -r bedfile; do
    outfile=$(basename "$bedfile" _peaks.narrowPeak)
    bedtools getfasta -fi "$REF_G38" -bed "$bedfile" -fo "fastas_GRCh38/${outfile}.fa"
    echo "$(realpath "fastas_GRCh38/${outfile}.fa")"
done < ref_GRCh38_bed_paths.txt > ref_GRCh38_fastas.txt

cat ref_GRCh37_fastas.txt ref_GRCh38_fastas.txt > all_fastas.txt
```

### Mode D analysis

#### Within-reference comparisons

```bash
hammock data/ref1_fastas.txt data/ref1_fastas.txt \
    --outprefix results/ref1_comparison \
    -k 20 -w 200 --precision 20

hammock data/ref2_fastas.txt data/ref2_fastas.txt \
    --outprefix results/ref2_comparison \
    -k 20 -w 200 --precision 20
```

#### Cross-reference comparison

```bash
hammock data/ref1_fastas.txt data/ref2_fastas.txt \
    --outprefix results/cross_reference \
    -k 20 -w 200 --precision 20

hammock data/all_fastas.txt data/all_fastas.txt \
    --outprefix results/combined \
    -k 20 -w 200 --precision 20
```

#### Parameter variations

```bash
for k in 15 20 25 30; do
    hammock data/all_fastas.txt data/all_fastas.txt \
        --outprefix results/combined_k${k} \
        -k $k -w 200 --precision 20
done

for w in 100 200 300 400; do
    hammock data/all_fastas.txt data/all_fastas.txt \
        --outprefix results/combined_w${w} \
        -k 20 -w $w --precision 20
done
```

### Expected results

1. **Within-reference consistency**: similarity of replicates on the same reference.
2. **Cross-reference concordance**: same biological sample across references.
3. **Reference-specific artifacts**: whether reference choice shifts similarity.

### Analysis notes

- Mode D uses k-mer content in sequences; extracted peak sequences still depend on alignment and reference.
- Different references can change peak locations, widths, and extracted sequence content.
