import gffpandas.gffpandas as gffpd
import os
import argparse

parser = argparse.ArgumentParser(description='Parse GFF file and create gene bed files')
parser.add_argument('gff_file', help='Path to the GFF file')
parser.add_argument('output_dir', help='Directory to output the gene bed files')

args = parser.parse_args()
gff_file = args.gff_file
output_dir = args.output_dir

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Read the GFF file into a DataFrame
annotation = gffpd.read_gff3(gff_file)  

df = annotation.df

# Create a copy of the filtered data to avoid modifying the original
gene_entries = df[df['type'] == 'gene'].copy()
gene_entries['ID'] = gene_entries['attributes'].str.extract(r'ID=([^;]+)')
gene_to_location = dict(zip(gene_entries['ID'], (gene_entries['seq_id'],gene_entries['start'], gene_entries['end'])))

transcript_entries = df[df['type'] == 'transcript'].copy()
transcript_entries['ID'] = transcript_entries['attributes'].str.extract(r'ID=([^;]+)')
transcript_to_location = dict(zip(transcript_entries['ID'], zip(transcript_entries['seq_id'], transcript_entries['start'], transcript_entries['end'])))
transcript_entries['parent'] = transcript_entries['attributes'].str.extract(r'Parent=([^;]+)')
# Create dictionary mapping transcripts to their parent genes
transcript_to_gene = dict(zip(transcript_entries['ID'], transcript_entries['parent']))
# Create dictionary mapping parent genes to their transcript IDs
gene_to_transcripts = {}
for transcript_id, gene_id in transcript_to_gene.items():
    if gene_id not in gene_to_transcripts:
        gene_to_transcripts[gene_id] = []
    gene_to_transcripts[gene_id].append(transcript_id)

# Exons
exon_entries = df[df['type'] == 'exon'].copy()
exon_entries['ID'] = exon_entries['attributes'].str.extract(r'ID=([^;]+)')
exon_to_location = dict(zip(exon_entries['ID'], zip(exon_entries['seq_id'],exon_entries['start'], exon_entries['end'])))

exon_entries['parent'] = exon_entries['attributes'].str.extract(r'Parent=([^;]+)')

# Create dictionary mapping exons to their parent transcripts
exon_to_transcript = dict(zip(exon_entries['ID'], exon_entries['parent']))
# Create dictionary mapping parent transcripts to their exon IDs
transcript_to_exons = {}
for exon_id, transcript_id in exon_to_transcript.items():
    if transcript_id not in transcript_to_exons:
        transcript_to_exons[transcript_id] = []
    transcript_to_exons[transcript_id].append(exon_id)


# Finally, write a bedfile for each gene containging the exon locations
for gene_id, transcripts in gene_to_transcripts.items():
    #open bedfile with name gene_id.bed
    with open(f"{output_dir}/{gene_id}.bed", "w") as f:
        for transcript_id in transcripts:
            exons = transcript_to_exons[transcript_id]
            for exon_id in exons:
                location = exon_to_location[exon_id]
                f.write(f"{location[0]}\t{location[1]}\t{location[2]}\t{exon_id}\n")
