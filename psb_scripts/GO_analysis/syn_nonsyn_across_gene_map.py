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
# PSB
####################################################################

from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



# below  are low counts codes for PSB 10 SNPs or more for 1 kb regions
gene_code = ['P0A8A4', 'A0P8X0', 'D0ZRB2', 'Q83ES6', 'P33941', 'P37877', 'Q81ST8', 'Q8ZNY9', 'D4G3R4', 'O34739', 'P46481', 'Q0VZ70', 'P76052', 'P76149', 'Q0P9C7', 'P23837', 'P45745', 'P69348', 'Q9KJ21', 'P9WN57', 'P76562', 'P75990', 'P0AES0', 'P9WMJ1', 'O87941', 'Q9KJ22', 'P0A734', 'P07061', 'P17054', 'P26647']
gene_codes = ['rhsC_4', 'wapA_3', 'rapA_5','rapA_12','rapA_9','P0AEZ3','P33919','C0SP86','P60240', 'hsdR_4', 'A0A0H2VDN9', 'D0ZRB2', 'E0SGL7', 'A1IGV8', 'Q9KIG4', 'G4NYJ6']
gene_names = [ 'Putative deoxyribonuclease RhsC', 'tRNA3(Ser)-specific nuclease WapA_3', 'RNA polymerase-associated protein RapA_5','RNA polymerase-associated protein RapA_12','RNA polymerase-associated protein RapA_9','Septum site-determining protein MinD', 'Putative DNA repair helicase RadD','DNA translocase SftA', 'RNA polymerase-associated protein RapA', 'Type-1 restriction enzyme R protein hsdR',  'Secretory immunoglobulin A-binding protein EsiB', 'E3 ubiquitin-protein ligase SlrP', 'Putative deoxyribonuclease RhsC', 'Alpha-agarase agaA33', 'Serine/threonine-protein kinase PK-1 spk1', 'tRNA3(Ser)-specific nuclease WapA']
gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# Example gff entry
# psb-scaff03    Prodigal:2.6    CDS    1660    2769    .    +    0    ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Get the information about start and end positions for each gene of interest
gene_location = {}
for n in range(0,len(gene_codes)):
    gene_code = gene_codes[n]
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
gene_name, seq = 'wapA_3', sequences['wapA_3']

# Reset the syn / nonsyn counters for each gene
nonsynonymous = 0
synonymous = 0
snp_counter = 0
protein_seq = ''

# Check which SNP positions are contained in the particular sequence
snp_positions = []
snp_frequency = []
for i in range(0, len(positions)):
    # if position of SNP is in the region
    if int(positions[i]) >= gene_location[gene_name][0] and int(positions[i]) <= gene_location[gene_name][1]:
        snp_positions.append(int(positions[i]))
        snp_counter += 1

        # Calculate SNP frequency: How many berries have this mutation?
        tmp = np.array(data[positions[i]])
        snp_frequency.append(np.nansum(tmp)/np.count_nonzero(~np.isnan(tmp)))

# Check if we need to reverse the strand. If so, reverse it.
if gene_location[gene_name][2] == '-':
    new_seq = ''
    for i in seq[::-1]:
        new_seq += bases[i]
    seq = new_seq

# Translate DNA to protein. Look for amino acid with a SNP
snp_positions = [i-gene_location[gene_name][0] for i in snp_positions]
nonsyn_position = []
syn_position = []
for i in range(0, len(seq), 3):     # For each codon

    # Get the original codon and amino acid
    protein = codon_aa[seq[i:i+3]]
    protein_seq += protein
    
    for j in range(0, len(snp_positions)):		# For each snp found in the region
        snp = snp_positions[j]
        #print(snp)
        if snp in range(i, i+3):	# If we're on the good codon

            # Get the mutant amino acid
            pos = snp % 3
            mutant = list(seq[i:i+3])
            mutant[pos] = mutations[snp+gene_location[gene_name][0]]
            mutant = ''.join(mutant)
            mutant_aa = codon_aa[mutant]

            #if snp < 600:
            #     print(protein + ' ' + mutant_aa + ' ' + str(int(snp/3)))

            # Categorize the change as synonymous / nonsynonymous
            if protein != mutant_aa:
                nonsynonymous += 1
                nonsyn_position.append(snp_frequency[j])  # add some jitter
                syn_position.append(np.nan)
                #print('nonsynonymous %d' %(snp+reg[0]))
            elif protein == mutant_aa:
                synonymous += 1
                nonsyn_position.append(np.nan)
                syn_position.append(snp_frequency[j])  # add some jitter
                #print('synonymous %d' %(snp+reg[0]))
            else:
                print('something went wrong')

print(protein_seq)
# Plot the mutations across the gene length
snp_positions = np.array(snp_positions)/3
plt.scatter(snp_positions, nonsyn_position, color='green', alpha=0.8, label='nonsynonymous')
plt.scatter(snp_positions, syn_position, color='red', alpha=0.8, label='synonymous')
plt.title(gene_names[gene_codes.index(gene_name)] + ' mutations across the gene')
plt.xlabel('gene length (in amino acids)')
plt.ylabel('SNP frequency')
plt.legend()

plt.show()



