
#!/bin/bash - 
#===================================/data/blangme2/jessica/remap2022/============================================
#
#          FILE: generate.sh
# 
#         USAGE: ./generate.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 05/06/2024 18:32
#      REVISION:  ---
#===============================================================================

# wget https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz

# wget https://remap.univ-amu.fr/storage/remap2022/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz

# gunzip remap2022_crm_macs2_hg38_v1_0.bed.gz remap2022_crm_macs2_hg19_v1_0.bed.gz

ln -s /data/blangme2/jessica/remap2022/remap2022_crm_macs2_hg38_v1_0.bed 
ln -s /data/blangme2/jessica/remap2022/remap2022_crm_macs2_hg19_v1_0.bed 
mkdir /data/blangme2/jessica/remap2022/subbeds
ln -s /data/blangme2/jessica/remap2022/subbeds

bedfile1="remap2022_crm_macs2_hg38_v1_0.bed"
# bedfile2="remap2022_crm_macs2_hg19_v1_0.bed"
proteins1=$(sed 's/\t/|/g' $bedfile1 | cut -d'|' -f4 | sort -u)
# proteins2=$(sed 's/\t/|/g' $bedfile2 | cut -d'|' -f4 | sort -u)
# proteins=$(echo $proteins1,$proteins2 | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
proteins=$(echo $proteins1 | tr ',' '\n' | tr ' ' '\n'| sort -u | tr '\n' ',' | sed 's/,$//')
echo $proteins | tr ',' '\n' > proteins.txt


for name in $(echo $proteins| tr ',' ' '); do
#   grep $name $bedfile1 > subbeds/${name}_hg38.bed
  awk -v name=$name -v OFS="\t" '{print $1,$2,$3}' subbeds/${name}_hg38.bed > subbeds/${name}_hg38_c3.bed
  
#   grep $name $bedfile2 > subbeds/${name}_hg19.bed
done
find subbeds -type f -empty -delete

# for name in {RAD21,SMC1A,SMC3,STAG1,STAG2,CTCF,YY1,ZEB1,STAT3,KMT2A,CREB1,GTF2F1,ESR1,REPIN1.ZNF157}; do
#   grep $name $bedfile1 > subbeds/${name}_hg38.bed
# done

find $PWD/subbeds/*hg38.bed -maxdepth 1 -type f > TFbeds.txt
find $PWD/subbeds/*hg38_c3.bed -maxdepth 1 -type f > TFbeds_c3.txt

# cd ../results
# grep -E "RAD21|SMC1A|SMC3|STAG1|STAG2|CTCF|YY1|ZEB1|STAT3|KMT2A|CREB1|GTF2F1|ESR1|REPIN1.ZNF157|TRIM22" ../data/TFbeds_c3.txt > ../data/TFbeds_c3_sub.txt
#  python3 ../../deterministic/lib/bed_minhash.py ../data/TFbeds_c3.txt 600 C
#  python3 ../../deterministic/lib/bed_jaccards_parallel.py ../data/TFbeds_c3_sub.txt C actual
#  head -n1 minhash_h600_matrixC.csv | sed 's/,/\n/g' | grep -E -n "RAD21|SMC1A|SMC3|STAG1|STAG2|CTCF|YY1|ZEB1|STAT3|KMT2A|CREB1|GTF2F1|ESR1|REPIN1.ZNF157|TRIM22"
# 111:CREB1_hg38_c3.bed
# 124:CTCFL_hg38_c3.bed
# 125:CTCF_hg38_c3.bed
# 177:ESR1_hg38_c3.bed
# 246:GTF2F1_hg38_c3.bed
# 358:KMT2A_hg38_c3.bed
# 575:RAD21_hg38_c3.bed
# 662:SMC1A-B_hg38_c3.bed
# 663:SMC1A_hg38_c3.bed
# 665:SMC3_hg38_c3.bed
# 708:STAG1_hg38_c3.bed
# 709:STAG2_hg38_c3.bed
# 712:STAT3_hg38_c3.bed
# 776:TRIM22_hg38_c3.bed
# 809:YY1AP1_hg38_c3.bed
# 810:YY1_hg38_c3.bed
# 840:ZEB1_hg38_c3.bed


 awk -v FS="," 'NR==1{print $1,$125,$575,$663,$665,$708,$776};NR==125{print $1,$125,$575,$663,$665,$708,$776};NR==575{print $1,$125,$575,$663,$665,$708,$776};NR==663{print $1,$125,$575,$663,$665,$708,$776};NR==665{print $1,$125,$575,$663,$665,$708,$776};NR==708{print $1,$125,$575,$663,$665,$708,$776};NR==776{print $1,$125,$575,$663,$665,$708,$776}' minhash_h600_matrixC.csv > CTCF_minhash_jaccardB.csv

awk -v FS="," 'NR==1{print $0};NR==125{print $0}' minhash_h600_matrixC.csv > CTCF_minhash_jaccard.csv

awk -v FS="," 'NR==1{print $0}' minhash_h600_matrixC.csv | sed 's/,/\n/g'> header_minhash_jaccard.txt

awk -v FS="," 'NR==125{print $0}' minhash_h600_matrixC.csv | sed 's/,/\n/g' > CTCF.tmp 
paste header_minhash_jaccard.txt CTCF.tmp > CTCF_minhash_jaccard.txt
sort -k2 -n -r CTCF_minhash_jaccard.txt > CTCF_minhash_jaccard_sorted.txt

