## Mus-Homo DNASE Experiments from Report Captured July 31, 2025

ENCODE dnase-Seq experiments were filtered for those tissues for which there were experiment results for BOTH Mus musculus and Homo sapiens. This yeilded a report (data/experiment_report_2025_7_31_19h_28m.tsv) and a text file fed to xargs to download a set of bedfiles each labeled by file accession number, found in the report under the Files column.


## Data Treatment
This code assumes you are in your data directory
### Download Data
```
tag=""
 mkdir -p beds${tag}
 cd beds${tag}
 xargs -n 1 curl -O -L < ../20250731_dnase_files.txt
cd ..
ls beds${tag}/* | xargs realpath > bed_paths${tag}.txt
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
ENCODE_report_to_key.sh experiment_report_2025_7_31_19h_28m.tsv accession_key_2025_7_31.tsv

# Filter it for what is actually there
   awk 'FNR==1 && FILENAME==ARGV[2] {print; next} FNR==NR {files[$0]=1; next} $2 in files' file_accessions.txt accession_key_2025_7_31.tsv > filtered_accession_key.tsv

# Create a more balanced balanced_accession_key.txt (balanced between species)
./balance_species.sh filtered_accession_key.tsv balanced_accession_key.tsv

# Filter for species
grep "Homo sapiens" filtered_accession_key.tsv | cut -f2 > Homo_sapiens_accessions.txt
grep -f Homo_sapiens_accessions.txt bed_paths.txt > Homo_sapiens_paths_beds.txt

grep "Homo sapiens" balanced_accession_key.tsv | cut -f2 > Homo_sapiens_accessions_balanced.txt
grep -f Homo_sapiens_accessions_balanced.txt bed_paths.txt > Homo_sapiens_paths_beds_balanced.txt


grep "Mus musculus" filtered_accession_key.tsv | cut -f2 > Mus_musculus_accessions.txt
grep -f Mus_musculus_accessions.txt bed_paths.txt > Mus_musculus_paths_beds.txt

grep "Mus musculus" balanced_accession_key.tsv | cut -f2 > Mus_musculus_accessions_balanced.txt
grep -f Mus_musculus_accessions_balanced.txt bed_paths.txt > Mus_musculus_paths_beds_balanced.txt
```

### Create Fasta Files using appropriate reference genomes

```
# Create directories for fasta files
mkdir -p fastas${tag}
# NOTE: instead  of making directory you could make a soft link

mm10=/data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta 

# Process mouse files using mm10 reference
while read -r bedfile; do
    # Remove both .bed.gz extensions
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $mm10 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa"
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Mus_musculus_paths_beds.txt > Mus_musculus_fastas.txt

grep -f Mus_musculus_accessions_balanced.txt Mus_musculus_fastas.txt > Mus_musculus_fastas_balanced.txt


grc38=/data/blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


# Process human files using hg38 reference  
while read -r bedfile; do
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $grc38 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa";
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Homo_sapiens_paths_beds.txt | sed 's/^[^/]*//' > Homo_sapiens_fastas.txt

grep -f Homo_sapiens_accessions_balanced.txt Homo_sapiens_fastas.txt | sed 's/^[^/]*//' > Homo_sapiens_fastas_balanced.txt

cat Mus_musculus_fastas.txt Homo_sapiens_fastas.txt > mus_homo_fastas.txt

cat Mus_musculus_fastas_balanced.txt Homo_sapiens_fastas_balanced.txt | sed 's/^[^/]*//' > mus_homo_fastas_balanced.txt

```



## Mode A/B/C
From within the TF3 directory (where you find this README)

```
hammock data/Mus_musculus_paths_beds.txt data/Mus_musculus_paths_beds.txt --mode C --precision 20 --outprefix results/mouseonly

hammock data/Homo_sapiens_paths_beds.txt data/Homo_sapiens_paths_beds.txt --mode C --precision 20 --outprefix results/humanonly

hammock data/Mus_musculus_paths_beds_balanced.txt data/Mus_musculus_paths_beds_balanced.txt --mode C --precision 20 --outprefix results/mouseonly_balanced

hammock data/Homo_sapiens_paths_beds_balanced.txt data/Homo_sapiens_paths_beds_balanced.txt --mode C --precision 20 --outprefix results/humanonly_balanced


```


## Mode D

```
hammock data/Mus_musculus_fastas.txt data/Mus_musculus_fastas.txt --outprefix results/mouseonly -k 20 -w 200 --precision 20

hammock data/Homo_sapiens_fastas.txt data/Homo_sapiens_fastas.txt --outprefix results/humanonly -k 20 -w 200 --precision 20

hammock data/mus_homo_fastas.txt data/mus_homo_fastas.txt --outprefix results/manmouse -k 20 -w 200 --precision 20

hammock data/mus_homo_fastas_balanced.txt data/mus_homo_fastas_balanced.txt --outprefix results/manmouse_balanced -k 10 -w 100 --precision 24

```

## Visualizations
Hammock 0.2.1 allows calling Rscript and includes a script to draw three heatmaps for a given output as long as a "report" can be provided in the same format as an ENCODE report. It takes the following arguments: `<ENCODE REPORT> <HAMMOCK OUTPUT> [<OUTPUT PREFIX>]` 

```{bash}
encode_heatmap.R data/experiment_report_2025_6_11_22h_5m.tsv results/manmouse_mnmzr_p20_jaccD_k20_w200.csv results/manmouse_p20k20w200jaccD


encode_heatmap.R data/experiment_report_2025_6_11_22h_5m.tsv results/manmouse_balanced_mnmzr_p24_jaccD_k10_w100.csv results/manmouse_balanced_p24k10w100jaccD

```

## Parameter Scan: Window, k, precision
First we seek out the best precision value when comparing mode B against bedtools jaccard as "ground truth," then we compare the results from bedtools jaccard on bedfiles against hammock. The process produces a lot of temp files, so we run it in a separate directory.

```
mkdir -p parameter_scan
cd parameter_scan
complete_parameter_sweep.sh -b ../data/bed_paths.txt -f ../data/mus_homo_fastas.txt -o scan




```