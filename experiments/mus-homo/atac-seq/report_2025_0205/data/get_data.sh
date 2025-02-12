#!/bin/bash
# experiment_report=experiment_report_2025_2_5_19h_35m.tsv
experiment_prefix=experiment_report_2025_2_5_20h_24m

# Find which columns have "Files"/"Accession" in the second row
files_column=$(awk -F'\t' 'NR==2 {for(i=1;i<=NF;i++) if($i=="Files") print i}' ${experiment_prefix}.tsv)

accession_column=$(awk -F'\t' 'NR==2 {for(i=1;i<=NF;i++) if($i=="Accession") print i}' ${experiment_prefix}.tsv)



experiment_report=${experiment_prefix}_no_header.tsv
tag="_n75"

# honestly that top row of the experiment report is a bit of a pain to parse
# so I'm going to remove it and the header row
awk -F'\t' 'NR>2' ${experiment_prefix}.tsv > ${experiment_report}

 mkdir -p processed${tag}
 cd processed${tag}
 xargs -n 1 curl -O -L < ../files_processed${tag}.txt
cd ..
ls processed${tag}/* | xargs realpath > processed_paths${tag}.txt


# For each experiment accession, get the list of files
mkdir -p experiment_files${tag}
# Store column 2 output in an array of experiment accessions
accession_list=($(awk -F'\t' -v col="$accession_column" '{print $col}' ${experiment_report}))
echo "${accession_list[@]}"

# Now you can iterate through the list of experiments to retrieve the lists of files
for thing in "${accession_list[@]}"; do
    grep -E "$thing" ${experiment_report} | awk -F'\t' -v col="$files_column" '{print $col}' | tr ',' '\n' | sed -E 's/^\/files\/(.*)\/$/\1/' > experiment_files${tag}/"$thing".txt
done

mkdir -p sublists${tag}
experiments=("Mus musculus" "Homo sapiens" "kidney" "liver" "lung" "cerebellum" "stomach")

for item in "${experiments[@]}"; do
    filename="${item// /_}"  # Replace spaces with underscores
    grep -E "$item" ${experiment_report} | awk -F'\t' -v col="$files_column" '{print $col}' | tr ',' '\n' | sed -E 's/^\/files\/(.*)\/$/\1/' | grep -v '^$' > sublists${tag}/"$filename"_files.txt

    grep -E "$item" ${experiment_report} | awk -F'\t' -v col="$accession_column" '{print $col}' | grep -v '^$' > sublists${tag}/"$filename"_experiments.txt

done


mkdir -p processed_paths${tag}

for item in "${experiments[@]}"; do
    rm -f processed_paths${tag}/"${item// /_}"_paths.txt
    filename=sublists${tag}/"${item// /_}"_files.txt
    # Get paths of processed files matching the pattern
    while read -r pattern; do
        
        find processed${tag} -type f -name "*$pattern*" -exec realpath {} \; >> processed_paths${tag}/"${item// /_}"_paths.txt
    done < "$filename"
done
    



#relink=https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz

