##############################################
# 
# dN/dS analysis SRB
#
##############################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def theory_dN(dS):
	 return (asymptotic_dNdS+(1-asymptotic_dNdS)*(1-np.exp(-sbymu*dS))/(theory_ds*sbymu))*dS

# Purifying selection model
# Fitted manually
asymptotic_dNdS = 0.12
dStar = 3e-04
sbymu = 1/dStar/asymptotic_dNdS
print("s/u =", sbymu)
print("s =", sbymu*1e-09)
theory_ds = np.logspace(-6,-1,100)


#theory_dNdSs = asymptotic_dNdS+(1-asymptotic_dNdS)*(1-numpy.exp(-sbymu*theory_ds))/(theory_ds*sbymu)
theory_dNdSs = theory_dN(theory_ds)/theory_ds

# Read the data.
data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_6_no_PB93.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()[1:]
pos = data.columns.tolist()
contigs = [str(i).split('.')[0] for i in data.iloc[0,:].tolist()]
data_positions = [pos[i].split('.')[0] + '_' + contigs[i] for i in range(0, len(pos))]   # make sure positions are whole numbers


orf_regions = set()
reg_to_descr = {}
coding_reg = [False] * len(data_positions)

print('finding all SNPs in ORFs.')
# Using the original GFF file, figure out which positions from the SNP matrix hit ORFs.
# Read in the annotation file.
with open('/home/ada/Desktop/Shraiman_lab/srb_data/otuB_mp_noG.gff', 'r') as f:
	for line in f:
		# Ignore the header.
		if not line.startswith('#') and not len(line.split('\t')) != 9: # Ignore header and sequence at the end of file.
			# We found a line with annotated gene. Check if any position from the matrix is inside the gene region.
			line = line.split('\t')
			contig, ann_start, ann_end, strand, descr = line[0].split('_')[1], int(line[3]), int(line[4]), line[6], line[8].split(';')[-1]
			for i in range(0, len(data_positions)):
				pos_tmp = int(data_positions[i].split('_')[0])
				contig_tmp = data_positions[i].split('_')[1]
				if pos_tmp >= ann_start and pos_tmp <= ann_end and contig_tmp == contig:
					# Position in the coding region AND in the correct contig -> we found an ORF.
					coding_reg[i] = True
					if (ann_end+1-ann_start)%3 == 0:
						# Make sure that we are looking at an actually transcribed protein and not an mRNA annotation.
						# (the length will be dividable by 3 (length of a codon))
						orf_regions.add((ann_start, ann_end, strand, contig))
						reg_to_descr[(ann_start, ann_end, strand, contig)] = descr
						break

print('Done!')
print('Fraction of coding positions where there are SNPs: %f' %(float(sum(coding_reg))/ len(coding_reg)))
sorted(orf_regions)

print('Create all ORF DNA sequences')
# Create a dictionary that will store the data about the ORF start/end positions and the actual FASTA sequence.
sequences = {}
with open('/home/ada/Desktop/Shraiman_lab/srb_data/srb_contigs.fasta', 'r') as f:
	current_contig = ''
	for line in f:
		if line.startswith('>'):
			current_contig = line.split('_')[1]
		else:
			for region in orf_regions:
				if region[3] == current_contig:
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
for i in data_positions:     # i is in "pos_contig" format
	mutations[i] = ''

# Check what was the mutation (minor allele) for each SNP position.
# Format: allele_positions, allele_contigs, minor allele
with open('/home/ada/Desktop/PinkBerry_scripts_paper/srb_scripts/dN_dS/srb_minor_alleles.txt', 'r') as f:
#with open('/home/ada/Desktop/PinkBerry_scripts_paper/srb_scripts/dN_dS/srb_minor_alleles_coverage_3.txt', 'r') as f:
	for line in f:
		pos, contig, min_allele = line.strip().split(',')
		tmp = pos + '_' + contig
		if tmp in mutations.keys():
			mutations[tmp] = min_allele


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
print('starting synonymous / nonsynonymous calculations for all possible loci')


positions_to_mut_type = {}
print(len(sequences.keys()))
counter = 0
for reg,seq in sequences.items():
	counter += 1
	print(counter)
	snp_positions = []
	for i in range(0, len(data_positions)):
		tmp_pos = int(data_positions[i].split('_')[0])
		tmp_contig = data_positions[i].split('_')[1]
		if tmp_pos >= int(reg[0]) and tmp_pos <= int(reg[1])and reg[3] == tmp_contig:
			snp_positions.append(tmp_pos)
	
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
				mutant[pos] = mutations[str(snp+reg[0])+'_'+reg[3]]
				mutant = ''.join(mutant)
				mutant_aa = codon_aa[mutant]

				# Categorize the change as synonymous / nonsynonymous
				if protein != mutant_aa:
					positions_to_mut_type[str(snp+reg[0])+'_'+reg[3]] = 'N'
				elif protein == mutant_aa:
					positions_to_mut_type[str(snp+reg[0])+'_'+reg[3]] = 'S'
				else:
					print('something went wrong')

print('Done!')		
#print([True for i in positions_to_mut_type.keys() if i in good_positions])

#data_index = data_index[:30]
#data_positions = data_positions[:2000]
print('Calculate experimental nonsynonymous and synonymous counts for each strain pair')
syn_differences = []
nonsyn_differences = []
all_syn_opportunities = []
all_nonsyn_opportunities = []
for m in range(0, len(data_index)):
	for n in range(m+1, len(data_index)):
		print(m, n)
		strain1 = data_index[m]
		strain2 = data_index[n]
		strain1_snps = np.array(data.loc[strain1, :])
		strain2_snps = np.array(data.loc[strain2, :])
		#print(strain1_snps)
		#print(strain2_snps)
		
		# Find differences between two strains 
		diff_SNPs = []
		for i in range(0, len(strain1_snps)):
			s1_snp = strain1_snps[i]
			s2_snp = strain2_snps[i]
			# Don't include sites where one of the positions has a NaN (we don't know what's there)
			if not np.isnan(s1_snp) and not np.isnan(s2_snp):
				if s1_snp != s2_snp:
					diff_SNPs.append(i)
		
		good_positions = [data_positions[i] for i in diff_SNPs]
		nonsynonymous, synonymous = 0, 0    # prepare a counter 
		for pos in good_positions:
			if pos in positions_to_mut_type.keys():
				if positions_to_mut_type[pos] == 'N':
					nonsynonymous += 1
				elif positions_to_mut_type[pos] == 'S':
					synonymous += 1

		# Find all genes that positions hit and count opportunity for synonymous and nonsynonymous changes.
		nonsyn_opportunity, syn_opportunity = 0, 0
		found_genes = []
		for reg,seq in sequences.items():
			for i in range(0, len(good_positions)):
				tmp_pos = int(good_positions[i].split('_')[0])
				tmp_contig = good_positions[i].split('_')[1]
				if tmp_pos >= int(reg[0]) and tmp_pos <= int(reg[1]) and tmp_contig == reg[3]:
					if reg not in found_genes:
						nonsyn_opportunity += gene_to_syn_nonsyn_opportunities[reg][1]   # value: [syn_count, nonsyn_count]
						syn_opportunity += gene_to_syn_nonsyn_opportunities[reg][0]
						found_genes.append(reg)

		syn_differences.append(synonymous)
		nonsyn_differences.append(nonsynonymous)
		all_syn_opportunities.append(syn_opportunity)
		all_nonsyn_opportunities.append(nonsyn_opportunity)

print('Done!')
print('Create a figure')
relative_dn = np.array(nonsyn_differences) / np.array(all_nonsyn_opportunities)
relative_ds = np.array(syn_differences) / np.array(all_syn_opportunities)
#experimental_dnds = np.array(nonsyn_differences) / np.array(syn_differences)
experimental_dnds = relative_dn / relative_ds

#plt.plot(theory_ds, theory_dNdSs, 'r-', label='Purifying selection model')
plt.plot(relative_ds, [1]*len(relative_ds), '--', color='gray', alpha=0.8, label='Neutral model')
plt.plot(relative_ds, experimental_dnds, 'o', color='0.7', markersize=2, markeredgewidth=0, label='(berry x berry)')
plt.xscale('log')
plt.yscale('log')
plt.title('Ratio of divergence at nonsynonymous sites ($d_N$) as a function of $d_S$  (SRB)')
plt.xlabel('Synonymous divergence, $d_S$')
plt.ylabel('Nonsynonymous ratio, $d_N/d_S$')
plt.show()
