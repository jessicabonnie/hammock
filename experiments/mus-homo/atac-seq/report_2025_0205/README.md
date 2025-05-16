## Mus-Homo ATAC-Seq Experiments from Report Captured Feb 25, 2025

ENCODE ATAC-Seq experiments were filtered for those tissues for which there were experiment results for BOTH Mus musculus and Homo sapiens. This yeilded a report (data/experiment_report_2025_2_5_20h_24m.tsv) and a text file fed to xargs to download a set of bedfiles each labeled by file accession number, found in the report under the Files column.

## Data Treatment
### Download Data
```
tag=""
 mkdir -p processed${tag}
 cd processed${tag}
 xargs -n 1 curl -O -L < ../files_processed${tag}.txt
cd ..
ls processed${tag}/* | xargs realpath > processed_paths${tag}.txt

```

### Make File Path Lists
```
tag=""
experiment_prefix=experiment_report_2025_2_5_20h_24m

# Find which columns have "Files"/"Accession" in the second row
files_column=$(awk -F'\t' 'NR==2 {for(i=1;i<=NF;i++) if($i=="Files") print i}' ${experiment_prefix}.tsv)

accession_column=$(awk -F'\t' 'NR==2 {for(i=1;i<=NF;i++) if($i=="Accession") print i}' ${experiment_prefix}.tsv)

experiment_report=${experiment_prefix}_no_header.tsv


awk -F'\t' 'NR>2' ${experiment_prefix}.tsv > ${experiment_report}

mkdir -p experiment_files${tag}
# Store column 2 output in an array of experiment accessions
accession_list=($(awk -F'\t' -v col="$accession_column" '{print $col}' ${experiment_report}))
echo "${accession_list[@]}"

# Now you can iterate through the list of experiments to retrieve the lists of files
for thing in "${accession_list[@]}"; do
    grep -E "$thing" ${experiment_report} | awk -F'\t' -v col="$files_column" '{print $col}' | tr ',' '\n' | sed -E 's/^\/files\/(.*)\/$/
\1/' > experiment_files${tag}/"$thing".txt
done

mkdir -p sublists${tag}
experiments=("Mus musculus" "Homo sapiens" "kidney" "liver" "lung" "cerebellum" "stomach")

for item in "${experiments[@]}"; do
    filename="${item// /_}"  # Replace spaces with underscores
    # Get the list of file accessions for each line in the report containing the keyword for the subset
    grep -E "$item" ${experiment_report} | awk -F'\t' -v col="$files_column" '{print $col}' | tr ',' '\n' | sed -E 's/^\/files\/(.*)\/$/\1/' | grep -v '^$' > sublists${tag}/"$filename"_files.txt

    # make a list of experiment accessions for each subgroup jic
    grep -E "$item" ${experiment_report} | awk -F'\t' -v col="$accession_column" '{print $col}' | grep -v '^$' > sublists${tag}/"$filename"_experiments.txt
done

# make lists of global file paths for subsets
mkdir -p processed_paths${tag}

for item in "${experiments[@]}"; do
    rm -f processed_paths${tag}/"${item// /_}"_paths.txt
    filename=sublists${tag}/"${item// /_}"_files.txt
    # Get paths of processed files matching the pattern
    while read -r pattern; do
        find processed${tag} -type f -name "*$pattern*" -exec realpath {} \; >> processed_paths${tag}/"${item// /_}"_paths.txt
    done < "$filename"
done

```

### Create Fasta Files using appropriate reference genomes

```
# Create directories for fasta files
mkdir -p data/fastas${tag}

mm10=/data/blangme2/jessica/mus_homo/references/mm10/mm10_no_alt_analysis_set_ENCODE.fasta 

# Process mouse files using mm10 reference
while read -r bedfile; do
    # Remove both .bed.gz extensions
    outfile=$(basename "$bedfile" .bed.gz)
    bedtools getfasta -fi $mm10 -bed "$bedfile" -fo "data/fastas${tag}/${outfile}.fa"
done < data/processed_paths${tag}/Mus_musculus_paths_beds.txt
cd data/fastas${tag}
realpath * > ../Mus_musculus_fastas.txt



# Process human files using hg38 reference  
while read -r bedfile; do
    outfile=$(basename "$bedfile" .bed)
    bedtools getfasta -fi /path/to/hg38/genome.fa -bed "$bedfile" -fo "fasta${tag}/${outfile}.fa"
done < processed_paths${tag}/Homo_sapiens_paths_beds.txt

# Create lists of fasta file paths
mkdir -p fasta_paths${tag}

<!-- # Create mouse fasta paths list
find fasta${tag} -type f -name "*$(cat processed_paths${tag}/Mus_musculus_paths.txt | xargs -I {} basename {} .bed)*.fa" \
    -exec realpath {} \; > fasta_paths${tag}/Mus_musculus_paths_fasta.txt -->

# Create human fasta paths list  
find fasta${tag} -type f -name "*$(cat processed_paths${tag}/Homo_sapiens_paths.txt | xargs -I {} basename {} .bed)*.fa" \
    -exec realpath {} \; > fasta_paths${tag}/Homo_sapiens_paths_fasta.txt




```



## Mode A/B/C
From within the report directory (where you find this README)

```
hammock data/Mus_musculus_paths_beds.txt data/Mus_musculus_paths_beds.txt --mode C --precision 20 --outprefix results/mouseonly
```


## Mode D

```
hammock data/Mus_musculus_fastas.txt data/Mus_musculus_fastas.txt --outprefix results/mouseonly

```