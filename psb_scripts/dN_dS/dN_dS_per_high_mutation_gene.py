####################################################################
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# For genes with high number of mutations, check their dN/dS value
# Is the dN/dS over 1.0 (positive selection)? Or is it similar to backgroun (less than 1.0 -> purifying selection)
# PSB
####################################################################


import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import pandas as pd


gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'
go_file_path = '/home/ada/Desktop/Shraiman_lab/GO_analysis/PSB_Uniprot_GO_all_genes.tsv'
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Example GO file entry
# P03004	P03004	DNAA_ECOLI	Chromosomal replication initiator protein DnaA	dnaA b3702 JW3679	Escherichia coli (strain K12)	DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]	ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]	cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; DnaA-DiaA complex [GO:1990102]; DnaA-Dps complex [GO:1990084]; DnaA-HU complex [GO:1990103]; DnaA-L2 complex [GO:1990082]; DnaA-oriC complex [GO:1990101]; plasma membrane [GO:0005886]; replication inhibiting complex [GO:1990078]; ATP binding [GO:0005524]; DNA binding [GO:0003677]; DNA replication origin binding [GO:0003688]; identical protein binding [GO:0042802]; sequence-specific DNA binding [GO:0043565]; DNA replication [GO:0006260]; DNA replication initiation [GO:0006270]; DNA unwinding involved in DNA replication [GO:0006268]; negative regulation of DNA-templated DNA replication initiation [GO:0032297]; positive regulation of DNA-templated DNA replication initiation [GO:0032298]; regulation of DNA replication [GO:0006275]; regulation of DNA-templated DNA replication initiation [GO:0030174]

# Example gff entry
# psb-scaff03	Prodigal:2.6	CDS	1660	2769	.	+	0	ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Read in the GO file and get all the gene codes
gene_codes = []
with open(go_file_path, 'r') as f:
	f.readline() # skip header
	for line in f:
		line = line.split('\t')
		gene_codes.append(line[1])


# Find all mutations in each gene
gene_length, gene_SNPs = [], []
gene_info = {}
counter = 0
for code in gene_codes:
	counter += 1
	# Specific to PSB data
	with open(gff_file_path, 'r') as f:
		for line in f:
			if line.startswith('psb-scaff03') and code in line:
				#print(counter)
				info = line.split('\t')[8]
				info = info.split(';')

				start, end = int(line.split('\t')[3]), int(line.split('\t')[4])

				gene_length.append(end-start)

				num_of_mutations = 0
				for i in data_positions:
					if int(i) >= start and int(i) <= end:
						num_of_mutations += 1
					elif int(i) > end:
						break

				gene_SNPs.append(num_of_mutations)
				key = str(start) + '_' + str(end)
				name = info[-1].split('=')[-1]
				value = [end-start, num_of_mutations, name]   # length, mutations, gene name
				gene_info[key] = value

# convert lists to np arrays
gene_length = np.array(gene_length)
gene_SNPs = np.array(gene_SNPs)

# Find a best fit line
a, b = np.polyfit(gene_length, gene_SNPs, 1)
std_SNPs = np.std(gene_SNPs)

# Check which genes have high number of mutations (higher than what we would expect from it based on length)
# Save those genes and check their dN/dS values for any genes which experience positive selection 

dnds_genes = []
for k,v in gene_info.items():
	# If the number of mutations is higher than the expected values + 1 standard deviation given the gene length
	if v[1] > (a*v[0]+b+std_SNPs):
		dnds_genes.append(k)

# Calculate dN/dS values for each gene
orf_regions = set()
print(len(dnds_genes))
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
					if (ann_end+1-ann_start)%3 == 0:
						# Make sure the gene is in the list of highly mutated genes to check
						if str(ann_start) + '_' + str(ann_end) in dnds_genes:  
							# Make sure that we are looking at an actually transcribed protein and not an mRNA annotation.
							# (the length will be dividable by 3 (length of a codon))
							orf_regions.add((ann_start, ann_end, strand))
							gene_order.append((ann_start, ann_end, strand))
							break

print('Done!')
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
print(len(data_positions))
with open('/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/dN_dS/psb_minor_alleles_coverage_3.txt', 'r') as f:
	for line in f:
		pos, minor = line.split(',')
		mutations[int(pos)] = minor.strip()

print('calculate opportunities: all possible synonymous and nonsynonymous changes for each gene')
gene_to_syn_nonsyn_opportunities = {}
good_positions = [int(i) for i in data_positions.copy()]
counter = 0
for reg,seq in sequences.items():
	counter += 1
	#print(counter)
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
#tmp_sequences = {list(sequences.keys())[8]: sequences[list(sequences.keys())[8]]}
#print(tmp_sequences)
#sequences = tmp_sequences
counter = 0
for reg,seq in sequences.items():
	counter += 1
	#print(counter)
	#print(reg)
	snp_positions = []
	for i in range(0, len(good_positions)):
		if int(good_positions[i]) >= int(reg[0]) and int(good_positions[i]) <= int(reg[1]):
			snp_positions.append(int(good_positions[i]))
	
	#print(snp_positions)
	# Check if we need to reverse the strand. If so, reverse it.
	if reg[2] == '-':
		new_seq = ''
		for i in seq[::-1]:
			new_seq += bases[i]
		seq = new_seq
	#print(seq)
	# Translate DNA to protein. Look for amino acid with a SNP
	snp_positions = [i-reg[0] for i in snp_positions]
	#print(snp_positions)
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
				#print(seq[i:i+3] + ' -> ' + protein + ' -> ' + mutations[snp+reg[0]] +' -> '+ mutant + ' -> ' + codon_aa[mutant])
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
counter = 0
genes_with_selection = []
dn_ds_values = []
gene_names_in_order = []

for i in range(0, len(gene_order)):
#for i in range(0, 600):
	print(str(i) + ' out of: ' + str(len(gene_order)))
	counter += 1
	current_gene = gene_order[i]
	
	# Find differences between two strains 
	diff_SNPs = []
	for j in range(0, len(data_positions)):
		if int(data_positions[j]) >= current_gene[0] and int(data_positions[j]) <= current_gene[1]:
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

	nonsyn_opportunity += gene_to_syn_nonsyn_opportunities[current_gene][1]   # value: [syn_count, nonsyn_count]
	syn_opportunity += gene_to_syn_nonsyn_opportunities[current_gene][0]

	syn_differences.append(synonymous)
	nonsyn_differences.append(nonsynonymous)
	all_syn_opportunities.append(syn_opportunity)
	all_nonsyn_opportunities.append(nonsyn_opportunity)

	relative_dn = nonsynonymous / nonsyn_opportunity
	relative_ds =  synonymous / syn_opportunity
	gene_length = current_gene[1] - current_gene[0]
	gene_names_in_order.append(gene_info[str(current_gene[0])+'_'+str(current_gene[1])][2].strip()+ str(counter) + '_' + str(len(good_positions)))
	if relative_dn / relative_ds >= 1.0:
		genes_with_selection.append(gene_info[str(current_gene[0])+'_'+str(current_gene[1])][2].strip()+ '_' + str(len(good_positions)))
		dn_ds_values.append(relative_dn / relative_ds)


print('Done!')
print('Create a figure')
relative_dn = np.array(nonsyn_differences) / np.array(all_nonsyn_opportunities)
relative_ds = np.array(syn_differences) / np.array(all_syn_opportunities)
experimental_dnds = relative_dn / relative_ds

#print(experimental_dnds)
#print(gene_names_in_order)

# Sort lists
sorted_indices = np.argsort(experimental_dnds)
gene_names_in_order = np.array(gene_names_in_order)
experimental_dnds_sorted = list(experimental_dnds[sorted_indices])
gene_names_in_order_sorted = list(gene_names_in_order[sorted_indices])

print(genes_with_selection)
print(dn_ds_values)

#print(experimental_dnds_sorted)
#print(gene_names_in_order_sorted)
#plt.barh(list(range(0,len(experimental_dnds_sorted))), experimental_dnds_sorted)
plt.barh(gene_names_in_order_sorted, experimental_dnds_sorted)
plt.axvline(x=1,color='red')
plt.title('dN/dS values for highly mutated genes')
plt.xlabel('dN/dS')
plt.show()
