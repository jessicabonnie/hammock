#!/bin/bash
# compare genes on chr1 and chr8 of diff reference should be close


mkdir gencode && cd gencode
# Obtain data fron gencode
#build 37
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.basic.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz

# gene_type protein coding and 3rd column "gene" --> 20,000 genes
zcat gencode.v46lift37.basic.annotation.gtf.gz | grep '^##' > filtered_lines_GRCh37.txt
zcat gencode.v46lift37.basic.annotation.gtf.gz | awk '$3 == "gene"' | grep 'gene_type "protein_coding"' >> filtered_genes_GRCh37.gtf

# build 38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# gene_type protein coding and 3rd column "gene" --> 20,000 genes
zcat gencode.v46.basic.annotation.gtf.gz |  grep '^##' >  filtered_genes_GRCh38.gtf
zcat gencode.v46.basic.annotation.gtf.gz | awk '$3 == "gene"' | grep 'gene_type "protein_coding"' >> filtered_genes_GRCh38.gtf

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz

# Example: Print "Hello, World!"
bedtools getfasta -fi GRCh37.primary_assembly.genome.fa -bed filtered_genes_GRCh37.gtf > gencode.v46lift37.BA.fasta 

bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed filtered_genes_GRCh38.gtf > gencode.v46.BA.fasta 

# End of your code