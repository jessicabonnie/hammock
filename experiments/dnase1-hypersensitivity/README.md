# DNase I Hypersensitivity Sites Analysis

## Overview

This experiment analyzes DNase I hypersensitivity sites (DHS) across different human fetal tissues to understand chromatin accessibility patterns and regulatory element conservation. The analysis focuses on computing pairwise similarity between tissue-specific DHS profiles using Jaccard indices and other similarity metrics.

## Background

DNase I hypersensitivity sites mark regions of open chromatin where regulatory proteins can access DNA. These sites are crucial indicators of active regulatory elements including promoters, enhancers, and other cis-regulatory sequences. 

### Dataset

The data used in this analysis comes from the Maurano et al. (2012) study published in Science: ["Systematic localization of common disease-associated variation in regulatory DNA"](https://www.science.org/doi/full/10.1126/science.1222794).

**Key characteristics of the dataset:**
- **Species**: Human fetal tissues  
- **Assay**: DNase I hypersensitivity sequencing (DNase-seq)
- **Processing**: Hotspot identification with two-pass algorithm, FDR < 0.05
- **Format**: BED files with merged hypersensitive regions
- **Tissues analyzed**: Brain, Heart, Small Intestine, Kidney, Lung, Muscle (arm/back/leg), Skin, Stomach
- **Total samples**: 20 fetal tissue samples across 8 tissue types

The dataset represents high-confidence DNase I hypersensitive sites identified using stringent statistical thresholds and merged to eliminate redundant overlapping regions.

## Methodology

### Data Processing Pipeline

1. **Input**: Pre-processed BED files containing merged DNase I hypersensitive sites
2. **Pairwise Analysis**: Computation of Jaccard similarity indices between all tissue pairs
3. **Similarity Matrix**: Generation of symmetric similarity matrix for downstream analysis

### Analysis Approach

The experiment uses hammock to compute pairwise Jaccard indices between all tissue samples, providing a quantitative measure of regulatory element overlap between different fetal tissues. It compares these results against those calculated by bedtools as seen in the [bedtools demo](http://quinlanlab.org/tutorials/bedtools.html)

## Code Used

### Data Acquisition

The dataset files can be obtained following the setup instructions from the [bedtools tutorial](http://quinlanlab.org/tutorials/bedtools.html):

```bash
# Navigate to the data directory
cd experiments/dnase1-hypersensitivity/data/

# Download the sample BED files
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/maurano.dnaseI.tgz
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/cpg.bed
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/exons.bed
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/gwas.bed
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/genome.txt
curl -O https://bedtools-tutorials.s3.amazonaws.com/curr-proc-bx/hesc.chromHmm.bed

# Extract the 20 DNase I hypersensitivity BED files
mkdir -p maurano.dnaseI
cd maurano.dnaseI/
tar -zxvf maurano.dnaseI.tgz
cd ..
rm maurano.dnaseI.tgz
# create an accession key from the filenames
echo -e echo -e "Accession\tFile\tTarget_of_Assay\tBiosample_term_name\tOrganism\tLife_stage" > maurano_filenames_key.tsv && ls maurano.dnaseI | awk -F'-' '{split($2,a,"\\."); print a[1] "\t" $0 "\t\t" $1 "\tHomo sapiens\tfetal"}' >> maurano_filenames_key.tsv

```

**Additional files included:**
- `cpg.bed`: CpG islands in the human genome
- `exons.bed`: RefSeq exons from human genes  
- `gwas.bed`: Human disease-associated SNPs from GWAS studies
- `hesc.chromHmm.bed`: Predicted chromatin states in human embryonic stem cells
- `genome.txt`: Genome coordinate file for analysis

These additional files enable comparative analyses between DNase I hypersensitive sites and other genomic features such as CpG islands, exonic regions, disease-associated variants, and chromatin states.

### Pairwise Jaccard Analysis

```bash
#!/bin/bash
# Compute pairwise Jaccard similarity for all DNase I hypersensitivity BED files

echo "file1 file2 intersection union jaccard n_intersections" > bedtools_pairwise_jaccard.txt

for file1 in `ls maurano.dnaseI/*.bed`
do
    for file2 in `ls maurano.dnaseI/*.bed`;
    do
        # Compute the jaccard statistic for these two files
        bvalues=`bedtools jaccard \
                   -a $file1 \
                   -b $file2 \
                   | awk 'NR==2 {print $0}'`

        echo $file1 $file2 $bvalues >> bedtools_pairwise_jaccard.txt
    done
    
done
```

## File Structure

```
experiments/dnase1-hypersensitivity/
├── README.md                     # This file
├── data/
│   ├── maurano.dnaseI/          # Source BED files (20 tissue samples)
│   │   ├── fBrain-DS14718.hotspot.twopass.fdr0.05.merge.bed
│   │   ├── fBrain-DS16302.hotspot.twopass.fdr0.05.merge.bed
│   │   ├── fHeart-DS15643.hotspot.twopass.fdr0.05.merge.bed
│   │   └── ... (17 more tissue BED files)
│   ├── cpg.bed                  # CpG islands in human genome
│   ├── exons.bed                # RefSeq exons from human genes
│   ├── gwas.bed                 # Human disease-associated SNPs
│   ├── hesc.chromHmm.bed        # Chromatin states in hESCs
│   ├── genome.txt               # Genome coordinates file
│   ├── maurano_files.txt        # List of input file paths
│   ├── maurano_files_key.tsv    # Key of accessions and biosamples etc.
│   ├── bedtools_pairwise.sh     # Pairwise analysis script
│   └── bedtools_pairwise_jaccard.txt  # Output similarity matrix
├── results/
│   └── maurano.R               # R analysis scripts
└── debug/                      # Debug and intermediate files
```

## Expected Results

The analysis generates a pairwise similarity matrix showing Jaccard indices between all tissue combinations. Key expected findings include:

- **Tissue-specific patterns**: Related tissues (e.g., different muscle samples) should show higher similarity
- **Developmental signatures**: Fetal tissues may show distinct accessibility patterns compared to adult tissues
- **Regulatory conservation**: Core regulatory elements may be conserved across multiple tissue types

## Usage

1. Set up the data directory and download files:
   ```bash
   # Create data directory if it doesn't exist
   mkdir -p experiments/dnase1-hypersensitivity/data/
   cd experiments/dnase1-hypersensitivity/data/
   
   # Follow data acquisition steps above
   ```

2. Run the pairwise analysis:
   ```bash
   # Follow pairwise jaccard analysis as described above
   ```

3. Analyze results:
   ```bash
   # View the similarity matrix
   head -10 ../bedtools_pairwise_jaccard.txt
   
   # Run R analysis (if available)
   Rscript ../../results/maurano.R
   
   # Optional: Use other tutorial files for comparative analysis
   bedtools intersect -a ../cpg.bed -b fBrain-DS14718.hotspot.twopass.fdr0.05.merge.bed
   ```

## Dependencies

- **bedtools**: Version 2.25+ for Jaccard analysis
- **R**: For downstream statistical analysis and visualization
- **Standard Unix tools**: awk, sed for text processing

## Notes

- All coordinates are in hg19/GRCh37 reference genome (inferred from `.hg19` in filenames)
- FDR threshold of 0.05 was used for hotspot identification (inferred from `.fdr0.05` in filenames)
- Merged regions eliminate overlapping DHS calls within samples (inferred from `.merge.bed` in filenames)
- Two-pass hotspot algorithm was used for processing (inferred from `.twopass` in filenames)

## Analysis

Create a fasta file for each bed file
```bash
grch37=/data/blangme2/fasta/grch37/hg19.fa
tag=""
# Create directories for fasta files
mkdir -p fastas${tag}
# Process human files using hg38 reference  
while read -r bedfile; do
    outfile=$(basename "$bedfile" .bed)
    bedtools getfasta -fi $grch37 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa";
    echo $(realpath fastas${tag}/${outfile}.fa)
done < maurano_files.txt > maurano_fastas.txt

```
Do a parameter sweep to look at how everything compares to the bedtools output for different parameter combos

```bash
# mkdir -p parameter_scan
# cd parameter_scan
# complete_parameter_sweep.sh -b ../data/maurano_files.txt -f ../data/maurano_fastas.txt -o maurano  -v --window 8,10,20,30,50,100,200 --klen 8,10,15,20,25 --precision 20,22,23,24 2> scan.log

#still in data folder
bfile=$(pwd)/maurano_files.txt
ffile=$(pwd)/maurano_fastas.txt

cd ..
mkdir -p cluster_analysis
cd cluster_analysis
ml parallel
sweepABC.sh $bfile

```


## Citation

If you use this analysis or data, please cite:

Maurano MT, Humbert R, Rynes E, et al. Systematic localization of common disease-associated variation in regulatory DNA. Science. 2012;337(6099):1190-1195. doi:10.1126/science.1222794
