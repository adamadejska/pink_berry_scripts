##############################################
# 
# dN/dS analysis PSB
# This script looks at a 10 gene windows and computes 
# the dN/dS for that window. At the end, it creates a 
# line plot showsing the changes in dN/dS based on genes
# across the chromosome. 
#
##############################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read the data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

orf_regions = set()
reg_to_descr = {}
coding_reg = [False] * len(data_positions)

print('finding all SNPs in ORFs.')
# Using the original GFF file, figure out which positions from the SNP matrix hit ORFs.
# Read in the annotation file.
gene_order = []
with open('/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff', 'r') as f:
	for line in f:
		# Ignore the header.
		if not line.startswith('#') and not len(line.split('\t')) != 9: # Ignore header and sequence at the end of file.
			# We found a line with annotated gene. Check if any position from the matrix is inside the gene region.
			line = line.split('\t')
			ann_start, ann_end, strand, descr = int(line[3]), int(line[4]), line[6], line[8].split(';')[-1]
			for i in range(0, len(data_positions)):
				if int(data_positions[i]) >= ann_start and int(data_positions[i]) <= ann_end:
					# Position in the coding region -> we found an ORF.
					coding_reg[i] = True
					if (ann_end+1-ann_start)%3 == 0:
						# Make sure that we are looking at an actually transcribed protein and not an mRNA annotation.
						# (the length will be dividable by 3 (length of a codon))
						orf_regions.add((ann_start, ann_end, strand))
						reg_to_descr[(ann_start, ann_end, strand)] = descr
						gene_order.append((ann_start, ann_end, strand))
						break

print('Done!')
print('Fraction of coding positions where there are SNPs: %f' %(float(sum(coding_reg))/ len(coding_reg)))
sorted(orf_regions)

print('Create all ORF DNA sequences')
# Create a dictionary that will store the data about the ORF start/end positions and the actual FASTA sequence.
sequences = {}
with open('/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_fasta_fixed.fasta', 'r') as f:
	for line in f:
		if not line.startswith('>') and not line.startswith('#'):
			for region in orf_regions:
				reg_s, reg_e = region[0], region[1]
				sequences[region] = line[reg_s-1:reg_e]
				
print('Done!')

# Dictionary with info about what codon codes for what amino acid. (* means STOP)
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
			'TAA':'*', 'TAG':'*', 'TGA':'*',
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

#  synonymous mutation is a change in the DNA sequence that codes for amino acids 
#  in a protein sequence, but does not change the encoded amino acid.
mutations = {}
for i in data_positions:
	mutations[int(i)] = ''

# Check what was the mutation (minor allele) for each SNP position.
with open('/home/ada/Desktop/Shraiman_lab/data/dsdn/minor_alleles.txt', 'r') as f:
	for line in f:
		tmp = line.split()
		for i in range(0, len(tmp)):
			mutations[int(data_positions[i])] = tmp[i][1:-1]

print('calculate opportunities: all possible synonymous and nonsynonymous changes for each gene')
gene_to_syn_nonsyn_opportunities = {}
good_positions = [int(i) for i in data_positions.copy()]
counter = 0
for reg,seq in sequences.items():
	counter += 1
	print(counter)
	nonsyn_count, syn_count = 0.0, 0.0
	for i in range(0, len(seq), 3):     # For each codon
		# Get the original codon and amino acid
		codon = seq[i:i+3]
		true_aa = codon_aa[codon]
		for p in range(0, 3):           # For each position in the codon
			# Find all possible changes at each position of the codon and how they will
			# affect the amino acid. 
			true_base = codon[p]
			nonsyn_tmp, syn_tmp = 0.0, 0.0
			for b in bases:
				if b != true_base:					
					tmp_aa = list(codon)
					tmp_aa[p] = b
					tmp_aa = codon_aa[''.join(tmp_aa)]
					if tmp_aa == true_aa:
						syn_tmp += 1
					else:
						nonsyn_tmp += 1
			
			nonsyn_count += nonsyn_tmp / 3
			syn_count += syn_tmp / 3

	gene_to_syn_nonsyn_opportunities[reg] = [syn_count, nonsyn_count]
			
print('Done!')
print('starting synonymous / nonsynonymous calculations for all possible experimental loci')

good_positions = [int(i) for i in data_positions.copy()]
positions_to_mut_type = {}
print(len(sequences.keys()))
counter = 0
for reg,seq in sequences.items():
	counter += 1
	print(counter)
	snp_positions = []
	for i in range(0, len(good_positions)):
		if int(good_positions[i]) >= int(reg[0]) and int(good_positions[i]) <= int(reg[1]):
			snp_positions.append(int(good_positions[i]))
	
	# Check if we need to reverse the strand. If so, reverse it.
	if reg[2] == '-':
		new_seq = ''
		for i in seq[::-1]:
			new_seq += bases[i]
		seq = new_seq

	# Translate DNA to protein. Look for amino acid with a SNP
	snp_positions = [i-reg[0] for i in snp_positions]
	for i in range(0, len(seq), 3):     # For each codon
		for snp in snp_positions:		# For each snp found in the region
			if snp in range(i, i+3):	# If we're on the good codon

				# Get the original codon and amino acid
				protein = codon_aa[seq[i:i+3]]

				# Get the mutant amino acid
				pos = snp % 3
				mutant = list(seq[i:i+3])
				mutant[pos] = mutations[snp+reg[0]]
				mutant = ''.join(mutant)
				mutant_aa = codon_aa[mutant]

				# Categorize the change as synonymous / nonsynonymous
				if protein != mutant_aa:
					positions_to_mut_type[snp+reg[0]] = 'N'
				elif protein == mutant_aa:
					positions_to_mut_type[snp+reg[0]] = 'S'
				else:
					print('something went wrong')

print('Done!')		

#data_positions = data_positions[:40]
print('Calculate experimental nonsynonymous and synonymous counts for each gene window')
syn_differences = []
nonsyn_differences = []
all_syn_opportunities = []
all_nonsyn_opportunities = []
gene_positions = []    # arbitrary
counter = 0
genes_with_selection = []
dn_ds_values = []

for i in range(0, len(gene_order)-10):
#for i in range(0, 600):
	print(str(i) + ' out of: ' + str(len(gene_order)))
	counter += 1
	current_genes = [gene_order[j] for j in range(i, i+10)]
	
	# Find differences between two strains 
	diff_SNPs = []
	for j in range(0, len(data_positions)):
		for gene in current_genes:
			if int(data_positions[j]) >= gene[0] and int(data_positions[j]) <= gene[1]:
				diff_SNPs.append(j)

	good_positions = [data_positions[j] for j in diff_SNPs]
	nonsynonymous, synonymous = 0, 0    # prepare a counter 
	for pos in good_positions:
		if int(pos) in positions_to_mut_type.keys():
			if positions_to_mut_type[int(pos)] == 'N':
				nonsynonymous += 1
			elif positions_to_mut_type[int(pos)] == 'S':
				synonymous += 1

	# Find all genes that positions hit and count opportunity for synonymous and nonsynonymous changes.
	nonsyn_opportunity, syn_opportunity = 0, 0
	for reg in current_genes:
		nonsyn_opportunity += gene_to_syn_nonsyn_opportunities[reg][1]   # value: [syn_count, nonsyn_count]
		syn_opportunity += gene_to_syn_nonsyn_opportunities[reg][0]

	syn_differences.append(synonymous)
	nonsyn_differences.append(nonsynonymous)
	all_syn_opportunities.append(syn_opportunity)
	all_nonsyn_opportunities.append(nonsyn_opportunity)
	gene_positions.append(counter)

	relative_dn = nonsynonymous / nonsyn_opportunity
	relative_ds =  synonymous / syn_opportunity
	if relative_dn / relative_ds >= 1.0:
		genes_with_selection.append(current_genes)
		dn_ds_values.append(relative_dn / relative_ds)


print('Done!')
print('Create a figure')
relative_dn = np.array(nonsyn_differences) / np.array(all_nonsyn_opportunities)
relative_ds = np.array(syn_differences) / np.array(all_syn_opportunities)
experimental_dnds = relative_dn / relative_ds

print(genes_with_selection)
print(dn_ds_values)

for i in range(0, len(experimental_dnds)-500, 500):
	tmp_experimental = experimental_dnds[i:i+500]
	tmp_neutral = [1]*len(relative_ds)
	tmp_neutral = tmp_neutral[i:i+500]
	plt.plot(tmp_neutral, '--', color='gray', alpha=0.8, label='Neutral model')
	plt.plot(tmp_experimental)
	#plt.title('Ratio of divergence at nonsynonymous sites ($d_N$) as a function of $d_S$')
	plt.xlabel('sliding 10-gene window')
	plt.ylabel('Nonsynonymous ratio, $d_N/d_S$')
	plt.show()
