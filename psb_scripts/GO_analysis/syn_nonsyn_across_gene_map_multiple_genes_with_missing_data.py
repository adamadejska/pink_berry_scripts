####################################################################
# We have a set of genes that are highly mutated in the PSB dataset.
# Are those mutations synonymous or nonsynonymous? 
# Nonsynonymous mutations might affect protein structure so we want
# to know how many are there and how they could affect the function. 
#
# How are those mutations distributed across the gene? 
# Each gene has different domains / repeats / regions with different functions.
# Are there nonsynonymous mutations in specific regions that would possibly
# change the functionality of the protein? 
#
# Some genes that code for toxins have immunity genes that are downstream.
# We want to see how the mutations look like in those regions as well.
# PSB
####################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


start, end = 4213680, 4234403

genes_wapA2 = [[3362633, 3366313], [3366315, 3367145], [3368130, 3368339], [3368740, 3370062], [3370053, 3370376], [3370712, 3372031], [3372010, 3373149], [3373426, 3375078]]
gene_names_wapA2 = ['wapA_2', 'hp1', 'hp2', 'hp3', 'hp4', 'hp5', 'hp6', 'hp7']

genes_wapA3 = [[5049505, 5058900], [5046919, 5049462], [5042095, 5046771], [5041002, 5041919], [5039083, 5041005], [5037116, 5039074]]
gene_names_wapA3 = ['wapA_3', 'hp1', 'hp2', 'hp3', 'hp4', 'hp5']

genes_wapA6 = [[6430294, 6435702], [6430132, 6430284], [6429077, 6429487], [6427834, 6428925], [6427067, 6427243]]
gene_names_wapA6 = ['wapA_6', 'hp1', 'hp2', 'hp3', 'hp4']

genes_rhsC3 = [[3350557, 3354828], [3354828, 3355520], [3355609, 3356232], [3356424, 3356621], [3356914, 3357006], [3357003, 3357290], [3357428,3357604], [3358062, 3358802], [3358855, 3359484]]
gene_names_rhsC3 = ['rhsC_3', 'hp1', 'hp2', 'hp3', 'hp4', 'hp5', 'hp6', 'hp7', 'hp8']

genes_rhsC4 = [[4219239, 4234403], [4218201, 4218800], [4217794, 4218252], [4217003, 4217176], [4214527, 4216869], [4214156, 4214494], [4213680, 4214045]]
gene_names_rhsC4 = ['rhsC_4', 'hp1', 'hp2', 'hp3', 'hp4', 'hp5', 'hp6']

gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'
current_genes = genes_rhsC4
current_names = gene_names_rhsC4
title = 'rhsC_4 and the hypothetical proteins proceeding it'

# Example gff entry
# psb-scaff03    Prodigal:2.6    CDS    1660    2769    .    +    0    ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Make a simple gene map without figuring out the kinds of mutations
# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
positions = data.columns.tolist()

snp_frequency = []
snp_positions = []
for i in positions:
    if int(i) >= start and int(i) <= end:
        # Calculate SNP frequency: How many berries have this mutation?
        tmp = np.array(data[i])
        snp_frequency.append(np.nansum(tmp)/np.count_nonzero(~np.isnan(tmp)))
        snp_positions.append(int(i))

color = ['gray'] * len(snp_positions)

# Get the information about start and end positions for each gene of interest
gene_location = {}
for n in range(0,len(current_genes)):
    gene_pos = current_genes[n]
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith('psb-scaff03'):
                    if str(gene_pos[0]) in line:
                        info = line.split('\t')[8]
                        info = info.split(';')

                        start, end, strand = int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6]

                        # Save the location for each gene
                        gene_location[current_names[n]] = [start, end, strand]



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

for gene_name in current_names:

    # Go through all the sequences to check for mutations
    seq = sequences[gene_name]

    # Reset the syn / nonsyn counters for each gene
    nonsynonymous = 0
    synonymous = 0
    protein_seq = ''

    # Check which SNP positions are contained in the particular sequence
    gene_snp_positions = []
    for i in range(0, len(positions)):
        # if position of SNP is in the region
        if int(positions[i]) >= gene_location[gene_name][0] and int(positions[i]) <= gene_location[gene_name][1]:
            gene_snp_positions.append(int(positions[i]))

    # Check if we need to reverse the strand. If so, reverse it.
    if gene_location[gene_name][2] == '-':
        new_seq = ''
        for i in seq[::-1]:
            new_seq += bases[i]
        seq = new_seq

    # Translate DNA to protein. Look for amino acid with a SNP
    gene_snp_positions = [i-gene_location[gene_name][0] for i in gene_snp_positions]
    nonsyn_position = []
    syn_position = []
    for i in range(0, len(seq), 3):     # For each codon

        # Get the original codon and amino acid
        protein = codon_aa[seq[i:i+3]]
        protein_seq += protein
        
        for j in range(0, len(gene_snp_positions)):		# For each snp found in the region
            snp = gene_snp_positions[j]

            if snp in range(i, i+3):	# If we're on the good codon

                # Get the mutant amino acid
                pos = snp % 3
                mutant = list(seq[i:i+3])
                mutant[pos] = mutations[snp+gene_location[gene_name][0]]
                mutant = ''.join(mutant)
                mutant_aa = codon_aa[mutant]

                snp_index = snp_positions.index(snp+gene_location[gene_name][0])
                # Categorize the change as synonymous / nonsynonymous
                if protein != mutant_aa:
                    nonsynonymous += 1
                    color[snp_index] = 'green'
                elif protein == mutant_aa:
                    synonymous += 1
                    color[snp_index] = 'red'
                else:
                    print('something went wrong')

# Divide the snp list into three lists so that we can have a legend
nonsyn = [np.nan] * len(snp_frequency)
syn = [np.nan] * len(snp_frequency)
noncoding = [np.nan] * len(snp_frequency)

for i in range(0, len(snp_frequency)):
    if color[i] == 'green':
        nonsyn[i] = snp_frequency[i]
    elif color[i] == 'red':
        syn[i] = snp_frequency[i]
    else:
         noncoding[i] = snp_frequency[i]

# calculate the fraction of missing data per locus
nan_frequency = []
for i in snp_positions:
    # Calculate NaN frequency: How much data is missing per locus?
    tmp = np.array(data[str(i)])
    nan_frequency.append(np.count_nonzero(np.isnan(tmp))/len(tmp))

plt.plot(snp_positions, nan_frequency, color='gray', label='missing data')

# Plot the mutations across the gene length
plt.scatter(snp_positions, nonsyn, color='green', alpha=0.8, label='nonsynonymous')
plt.scatter(snp_positions, syn, color='red', alpha=0.8, label='synonymous')
plt.scatter(snp_positions, noncoding, color='gray', alpha=0.8, label='noncoding')

# Color the area of each ORF on the graph
for orf in current_genes:
    plt.axvspan(orf[0], orf[1], facecolor='gray', alpha=0.3)
    plt.axvline(orf[0], color = 'gray', alpha=0.4)
    plt.axvline(orf[1], color = 'gray', alpha=0.4)

plt.xticks(range(snp_positions[0], snp_positions[-1], 1000))
plt.xlabel('Chromosomal position')
plt.ylabel('Frequency')
plt.title(title)
plt.legend()
plt.show()
