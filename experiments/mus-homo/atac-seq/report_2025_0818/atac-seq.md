## Mus-Homo DNASE Experiments from Report Captured July 31, 2025

ENCODE dnase-Seq experiments were filtered for those tissues for which there were experiment results for BOTH Mus musculus and Homo sapiens. This yeilded a report (data/experiment_report_2025_8_18_17h_30m.tsv) and a text file fed to xargs to download a set of bedfiles each labeled by file accession number, found in the report under the Files column.


## Data Treatment
This code assumes you are in your data directory
### Download Data
```
 mkdir -p beds
 cd beds
 xargs -n 1 curl -O -L < ../20250818_atacseq_files.txt
cd ..
ls beds/* | xargs realpath > bed_paths.txt
# make a list of the file accessions
# This handles .gz, .bz2, .xz and other extensions
while read -r filepath; do basename "$filepath" | sed -E 's/\.(gz|bz2|xz)$//' | sed 's/\.[^.]*$//'; done < bed_paths.txt > file_accessions.txt


# HUMAN
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.p14.genome.fa.gz
#lab location: /data/blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# MOUSE
#lab location: /data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta 

```

### Make A Key and Subset files as needed 
```
# Make a key
ENCODE_report_to_key.sh experiment_report_2025_8_18_17h_30m.tsv accession_key_2025_8_18.tsv
#sed 's/^[^/]*//' accession_key_2025_08_18.tsv > accession_key_2025_7_31_with_lifestage.tsv

# Filter it for what is actually there
   awk 'FNR==1 && FILENAME==ARGV[2] {print; next} FNR==NR {files[$0]=1; next} $2 in files' file_accessions.txt accession_key_2025_8_18.tsv > filtered_accession_key.tsv

# Create a table of tissues with sample counts
python summarize_accession_key.py filtered_accession_key.tsv
```
## All Available Samples

### Homo sapiens

**Total samples:** 6
**Unique tissues:** 5
**Life stages:** adult
**Life stage distribution:** adult (6, 100.0%)

| Tissue | Total Count | Percentage | Life Stage Breakdown |
|--------|-------------|------------|----------------------|
| stomach | 2 | 33.3% | adult: 2 |
| liver | 1 | 16.7% | adult: 1 |
| lung | 1 | 16.7% | adult: 1 |
| cerebellum | 1 | 16.7% | adult: 1 |
| kidney | 1 | 16.7% | adult: 1 |

### Mus musculus

**Total samples:** 19
**Unique tissues:** 5
**Life stages:** adult, embryonic, postnatal
**Life stage distribution:** adult (1, 5.3%), embryonic (14, 73.7%), postnatal (4, 21.1%)

| Tissue | Total Count | Percentage | Life Stage Breakdown |
|--------|-------------|------------|----------------------|
| liver | 7 | 36.8% | embryonic: 6, postnatal: 1 |
| stomach | 4 | 21.1% | embryonic: 3, postnatal: 1 |
| kidney | 4 | 21.1% | embryonic: 3, postnatal: 1 |
| lung | 3 | 15.8% | embryonic: 2, postnatal: 1 |
| cerebellum | 1 | 5.3% | adult: 1 |

```
# Filter for species
grep "Homo sapiens" filtered_accession_key.tsv | cut -f2 > Homo_sapiens_accessions.txt
grep -f Homo_sapiens_accessions.txt bed_paths.txt > Homo_sapiens_paths_beds.txt

grep "Mus musculus" filtered_accession_key.tsv | cut -f2  > Mus_musculus_accessions.txt
grep -f Mus_musculus_accessions.txt bed_paths.txt > Mus_musculus_paths_beds.txt

```

### Create Fasta Files using appropriate reference genomes

```
# Create directories for fasta files
mkdir -p fastas
# NOTE: instead  of making directory you could make a soft link to elsewhere

mm10=/data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta 

# Process mouse files using mm10 reference
while read -r bedfile; do
    # Remove both .bed.gz extensions
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $mm10 -bed "$bedfile" -fo "fastas/${outfile}.fa"
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Mus_musculus_paths_beds.txt > Mus_musculus_fastas.txt

sed -i 's/^[^/]*//' Mus_musculus_fastas.txt


## Human fastas

grc38=/data/blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


# Process human files using hg38 reference  
while read -r bedfile; do
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $grc38 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa";
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Homo_sapiens_paths_beds.txt > Homo_sapiens_fastas.txt


cat Mus_musculus_fastas.txt Homo_sapiens_fastas.txt > mus_homo_fastas.txt
```


## Subsets
Hammock includes a script to filter hammock output files to contain only the pairwise comparisons for a list of file accessions so that hammock does not need to be run again to have hammock output for only that subset.

```bash
mkdir mouse_
batch_filter_hammock.sh ../mouseonly/ ../data/embryonic_accessions.txt

```



## Visualizations
Hammock 0.2.1 allows calling Rscript and includes a script to draw three heatmaps for a given output as long as a "report" can be provided in the same format as an ENCODE report. It takes the following arguments: `<ENCODE REPORT> <HAMMOCK OUTPUT> [<OUTPUT PREFIX>]` 

```{bash}

encode_heatmap.R data/experiment_report_2025_7_31_19h_28m.tsv results/manmouse_balanced_mnmzr_p24_jaccD_k10_w100.csv results/manmouse_balanced_p24k10w100jaccD

encode_heatmap.R data/experiment_report_2025_7_31_19h_28m.tsv results/manmouse_subtissue_mnmzr_p24_jaccD_k10_w100.csv results/manmouse_subtissue_p24k10w100jaccD

```

## Parameter Scans: expA, subB, Window, k, precision

hammock 0.4.0 has a shiny app under `graphics/apps/cluster_analysis` that enables juxtaposing graphs from outputs from a scan of parameters.In order to create the inputs there are bash scripts in the scripts directory that can be called by name

```bash
mkdir -p mouseonly
cd mouseonly
sweepABC.sh  ../data/Mus_musculus_paths_beds.txt

sweepD.sh ../data/Mus_musculus_fastas.txt

```
The code to carry out the same analysis on single output files is included in `scripts/clustering_analysis.R`, `scripts/clustering_analysis.py`, which relies on other dependencies in that directory.

## Interpolate between A and B
This script will default to looking for precision 24 outputs in the provided directory. It will use whatever mode C outputs which have subB or expA parameters as well as the basic modeA run and modeB run with default names. It's lazy that way.

```bash
    # I'm assuming you were still in mouseonly directory ... but you know how to use relative paths
    cd ..
    interpolateABC.py mouseonly/ -t mouseonly 

```