import gffpandas.gffpandas as gffpd
import os
import argparse

parser = argparse.ArgumentParser(description='Parse GTF file and create gene bed files')
parser.add_argument('gtf_file', help='Path to the GTF file')
parser.add_argument('output_dir', help='Directory to output the gene bed files')

args = parser.parse_args()
gtf_file = args.gtf_file
output_dir = args.output_dir

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read the GFF file into a DataFrame
annotation = gffpd.read_gff3(gtf_file)  

df = annotation.df

# Create a copy of the filtered data to avoid modifying the original
gene_entries = df[df['type'] == 'gene'].copy()
gene_entries['ID'] = gene_entries['attributes']
gene_to_location = dict(zip(gene_entries['ID'], (gene_entries['seq_id'],gene_entries['start'], gene_entries['end'])))


# Exons
exon_entries = df[df['type'] == 'exon'].copy()
exon_entries['gene_id'] = exon_entries['attributes'].str.extract(r'gene_id "([^";]+)')

# Create a dictionary mapping gene IDs to a list of their exon position tuples
gene_to_exons = {}
for index, row in exon_entries.iterrows():
    gene_id = row['gene_id']
    if gene_id not in gene_to_exons:
        gene_to_exons[gene_id] = []
    gene_to_exons[gene_id].append((row['seq_id'], row['start'], row['end']))

# Write the exon locations to a bed file for each gene
for gene_id, exon_locations in gene_to_exons.items():
    with open(f"{output_dir}/{gene_id}.bed", "w") as f:
        for location in exon_locations:
            f.write(f"{location[0]}\t{location[1]}\t{location[2]}\n")