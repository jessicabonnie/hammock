## Mus-Homo TF3 Binding Experiments from Report Captured June 11, 2025

ENCODE Chip-Seq experiments were filtered for those tissues and Target TFs for which there were experiment results for BOTH Mus musculus and Homo sapiens. This yeilded a report (data/experiment_report_2025_6_11_22h_5m.tsv) and a text file fed to xargs to download a set of bedfiles each labeled by file accession number, found in the report under the Files column.

The three Transcription Factors with experiments in both species are CTCF (also analyzed separately), POLR2A, EP300. There were

## Data Treatment
This code assumes you are in your data directory
### Download Data
```
tag=""
 mkdir -p beds${tag}
 cd beds${tag}
 xargs -n 1 curl -O -L < ../20250611_TF3_files.txt
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
../../../scripts/ENCODE_report_to_key.sh experiment_report_2025_6_11_22h_5m.tsv accession_key_2025_6_11.tsv

# Filter it for what is actually there
awk 'NR==1 {print; next} FNR==NR {files[$0]=1; next} $2 in files' file_accessions.txt accession_key_2025_6_11.tsv > filtered_accession_key.tsv

# Filter for species
grep "Homo sapiens" filtered_accession_key.tsv | cut -f2 > Homo_sapiens_accessions.txt
grep -f Homo_sapiens_accessions.txt bed_paths.txt > Homo_sapiens_paths_beds.txt

grep "Mus musculus" filtered_accession_key.tsv | cut -f2 > Mus_musculus_accessions.txt
grep -f Mus_musculus_accessions.txt bed_paths.txt > Mus_musculus_paths_beds.txt

# Filter for TF Target

grep "POLR2A" filtered_accession_key.tsv | cut -f2 > POLR2A_accessions.txt
grep -f POLR2A_accessions.txt bed_paths.txt > POLR2A_paths_beds.txt

grep "EP300" filtered_accession_key.tsv | cut -f2 > EP300_accessions.txt
grep -f EP300_accessions.txt bed_paths.txt > EP300_paths_beds.txt

grep "CTCF" filtered_accession_key.tsv | cut -f2 > CTCF_accessions.txt
grep -f CTCF_accessions.txt bed_paths.txt > CTCF_paths_beds.txt

```

### Create Fasta Files using appropriate reference genomes

```
# Create directories for fasta files
mkdir -p fastas${tag}

mm10=/data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta 

# Process mouse files using mm10 reference
while read -r bedfile; do
    # Remove both .bed.gz extensions
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $mm10 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa"
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Mus_musculus_paths_beds.txt > Mus_musculus_fastas.txt


grc38=/data/blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


# Process human files using hg38 reference  
while read -r bedfile; do
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $grc38 -bed "$bedfile" -fo "fastas${tag}/${outfile}.fa";
    echo $(realpath fastas${tag}/${outfile}.fa)
done < Homo_sapiens_paths_beds.txt > Homo_sapiens_fastas.txt

cat Mus_musculus_fastas.txt Homo_sapiens_fastas.txt > mus_homo_fastas.txt

```



## Mode A/B/C
From within the TF3 directory (where you find this README)

```
hammock data/Mus_musculus_paths_beds.txt data/Mus_musculus_paths_beds.txt --mode C --precision 20 --outprefix results/mouseonly

hammock data/Homo_sapiens_paths_beds.txt data/Homo_sapiens_paths_beds.txt --mode C --precision 20 --outprefix results/humanonly
```


## Mode D

```
hammock data/Mus_musculus_fastas.txt data/Mus_musculus_fastas.txt --outprefix results/mouseonly -k 20 -w 200 --precision 20

hammock data/Homo_sapiens_fastas.txt data/Homo_sapiens_fastas.txt --outprefix results/humanonly -k 20 -w 200 --precision 20

hammock data/mus_homo_fastas.txt data/mus_homo_fastas.txt --outprefix results/manmouse -k 20 -w 200 --precision 20
```

## Visualizations
Hammock 0.2.1 allows calling Rscript and includes a script to draw three heatmaps for a given output as long as a "report" can be provided in the same format as an ENCODE report. It takes the following arguments: `<ENCODE REPORT> <HAMMOCK OUTPUT> [<OUTPUT PREFIX>]` 

```{bash}
encode_heatmap.R data/experiment_report_2025_6_11_22h_5m.tsv results/manmouse_mnmzr_p20_jaccD_k20_w200.csv results/manmouse_p20jaccD

```