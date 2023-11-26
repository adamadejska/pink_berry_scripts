####################################################################
# Translate the nucleotide sequence of given genes into amino acids
#
# PSB
####################################################################

from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


genes = ['rhsC_4', 'wapA_3', 'rsgI3_1', 'lgrB', 'cbhB', 'sdrF', 'spk1', 'wapA_2', 'mscM_2', 'rhsC_3', 'wapA_6', 'rapA_9', 'esiB_3', 'afsR', 'rapA_1', 'hsdR_4', 'slrP', 'egtB_3', 'rapA_5', 'xerD_15', 'sftA', 'hsdM', 'xerD_16', 'pksN_2', 'pksN_1']
gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# Example gff entry
# psb-scaff03    Prodigal:2.6    CDS    1660    2769    .    +    0    ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Get the information about start and end positions for each gene of interest
gene_location = {}
for n in range(0,len(genes)):
    gene_code = genes[n]
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith('psb-scaff03'):
                    if gene_code in line:
                        info = line.split('\t')[8]
                        info = info.split(';')

                        start, end, strand = int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6]

                        # Save the location for each gene
                        gene_location[gene_code] = [start, end, strand]

# Create a dictionary that will store the data about the ORF start/end positions and the actual FASTA sequence.
sequences = {}
with open('/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_fasta_fixed.fasta', 'r') as f:
    for line in f:
        if not line.startswith('>') or not line.startswith("#"):  # Ignore the header of the first line of file
            for gene,location in gene_location.items():
                sequences[gene] = line[location[0]-1:location[1]]   # Get the start and end of the gene by slicing
                

# Dictionary with info about what codon codes for what amino acid.
codon_aa = {'TTT': 'F', 'TTC': 'F',
			'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
			'ATT':'I', 'ATC':'I', 'ATA':'I',
			'ATG':'M',
			'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
			'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
			'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
			'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
			'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
			'TAT':'Y', 'TAC':'Y',
			'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP',
			'CAT':'H', 'CAC':'H',
			'CAA':'Q', 'CAG':'Q',
			'AAT':'N', 'AAC':'N',
			'AAA':'K', 'AAG':'K',
			'GAT':'D', 'GAC':'D',
			'GAA':'E', 'GAG':'E',
			'TGT':'C', 'TGC':'C',
			'TGG':'W',
			'CGT':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
			'AGT':'S', 'AGC':'S',
			'GGT':'G', 'GGA':'G', 'GGC':'G', 'GGG':'G'}

bases = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}  # base pairs used for converting + to - strands.

# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
positions = data.columns.tolist()

#  synonymous mutation is a change in the DNA sequence that codes for amino acids 
#  in a protein sequence, but does not change the encoded amino acid.
mutations = {}
for i in positions:
	mutations[int(i)] = ''

# Check what was the mutation (minor allele) for each SNP position.
with open('/home/ada/Desktop/Shraiman_lab/data/dsdn/minor_alleles.txt', 'r') as f:
	for line in f:
		tmp = line.split()
		for i in range(0, len(tmp)):
			mutations[int(positions[i])] = tmp[i][1:-1]

# Go through all the sequences to check for mutations
out = open('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/high_mutations_genes_amino_acid_seq.fasta', 'w')
for gene_name in genes:
    seq = sequences[gene_name]

    protein_seq = ''
    snp_counter = 0

    # Check which SNP positions are contained in the particular sequence
    for i in range(0, len(positions)):
        # if position of SNP is in the region
        if int(positions[i]) >= gene_location[gene_name][0] and int(positions[i]) <= gene_location[gene_name][1]:
            snp_counter += 1

    # Check if we need to reverse the strand. If so, reverse it.
    if gene_location[gene_name][2] == '-':
        new_seq = ''
        for i in seq[::-1]:
            new_seq += bases[i]
        seq = new_seq


    for i in range(0, len(seq), 3):     # For each codon

        # Get the original codon and amino acid
        protein = codon_aa[seq[i:i+3]]
        protein_seq += protein

    out.write('>' + gene_name + ' number of mutations: ' + str(snp_counter) + '\n')
    out.write(protein_seq + '\n')

