# A program to get the pairwise jaccard distances between two lists of bed files

# Usage: ./bedtools.sh <file1> <file2>

# Get the pairwise jaccard distance between two lists of bed files
ml bedtools
# read the list of bed file paths from the first argument
file1=$1
# read each line of the first file into an array
file1_lines=()
while IFS= read -r line; do
    file1_lines+=("$line")
done < "$file1"

# read the list of bed file paths from the second argument
file2=$2
# read each line of the second file into an array
file2_lines=()
while IFS= read -r line; do
    file2_lines+=("$line")
done < "$file2"

# get the pairwise jaccard distance between the two lists of bed files
for file1_line in "${file1_lines[@]}"; do
    sort -k1,1 -k2,2n $file1_line > $file1_line.sorted.bed
    for file2_line in "${file2_lines[@]}"; do
        # get the jaccard distance between the two bed files
        if [ ! -f $file2_line.sorted.bed ]; then
            sort -k1,1 -k2,2n $file2_line > $file2_line.sorted.bed
        fi
        jaccard_distance=$(bedtools jaccard -a $file1_line.sorted.bed -b $file2_line.sorted.bed | tail -n 1 | cut -f 3)
        # Suppress output
        # echo "$file1_line $file2_line $jaccard_distance"
    done
    rm $file1_line.sorted.bed
done
for file2_line in "${file2_lines[@]}"; do
    rm $file2_line.sorted.bed
done




