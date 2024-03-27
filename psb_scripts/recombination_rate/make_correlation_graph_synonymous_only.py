#############################################################################
# Create a graph similar to https://www.nature.com/articles/s41592-018-0293-7
# Where we choose a length l and calculate the probability that we find a mutation
# at distance i+l given that there was a mutation at i
#############################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read the data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

###############################################################################
# First, figure out which loci create synonymous mutations.
# We need to filter all loci so that only synonymous are left and use those
# for the correlation calculations
###############################################################################

orf_regions = set()
reg_to_descr = {}
coding_reg = [False] * len(data_positions)
coding_snps = []

print('finding all SNPs in ORFs.')
# Using the original GFF file, figure out which positions from the SNP matrix hit ORFs.
# Read in the annotation file.
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
					coding_snps.append(data_positions[i])
					coding_reg[i] = True
					if (ann_end+1-ann_start)%3 == 0:
						# Make sure that we are looking at an actually transcribed protein and not an mRNA annotation.
						# (the length will be dividable by 3 (length of a codon))
						orf_regions.add((ann_start, ann_end, strand))
						reg_to_descr[(ann_start, ann_end, strand)] = descr
						break

print('Done!')
coding_SNPs_count = sum(coding_reg)
noncoding_SNPs = sum(np.array(coding_reg) == False)
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
			
print('Done!')
print('starting synonymous / nonsynonymous calculations for all possible loci')
coding_snps = list(set(coding_snps))
sorted(coding_snps)
good_positions = [int(i) for i in coding_snps]
synonymous_loci = []
nonsynonymous_loci = []
print(len(sequences.keys()))

counter = 0
snp_counter = 0
for reg,seq in sequences.items():
	counter += 1
	print(counter)
	snp_positions = []
	for i in range(0, len(good_positions)):
		if int(good_positions[i]) >= int(reg[0]) and int(good_positions[i]) <= int(reg[1]):
			snp_positions.append(int(good_positions[i]))
			snp_counter += 1
	
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
				if protein == mutant_aa:
					synonymous_loci.append(snp+reg[0])
				if protein != mutant_aa:
					nonsynonymous_loci.append(snp+reg[0])

print('Done!')
"""
print('coding SNPs: ' + str(coding_SNPs_count))
print('noncoding SNPs: ' + str(noncoding_SNPs))
print('synonymous SNPs: ' + str(len(synonymous_loci)))
print('nonsynonymous SNPs: ' + str(len(nonsynonymous_loci)))

yvals = [noncoding_SNPs, coding_SNPs_count, len(synonymous_loci), len(nonsynonymous_loci)]
xticks = ['noncoding', 'coding', 'synonymous', 'nonsynonymous']
plt.bar(range(0, len(yvals)), yvals)
plt.xticks(ticks=range(0,len(yvals)), labels=xticks)
plt.yscale('log')
plt.show()
"""
synonymous_loci = list(set(synonymous_loci))
synonymous_loci.sort()
synonymous_loci = [str(i) for i in synonymous_loci]
# Subset the data so that it only includes loci with synonymous mutations
data = data.loc[:, synonymous_loci]
data_index = data.index.values.tolist()
data_positions = synonymous_loci

all_sigmas = []
# Make a matrix of vectors. Each vector is a binary vector where 1 = different alleles at pos i, and 0 = the same alleles at pos i
for s1i in range(0, len(data_index)):
	for s2i in range(s1i+1, len(data_index)):
		print(s1i, s2i)
		s1, s2 = data_index[s1i], data_index[s2i]
		v1, v2 = np.array(data.loc[s1, :]), np.array(data.loc[s2, :])

		sigma = [s1+'_'+s2]
		for i in range(0, len(v1)):
			if v1[i] not in [0, 1] or v2[i] not in [0, 1]:   # not enough info for calculations
				sigma.append(np.nan)
			elif v1[i] == v2[i]:                     # both alleles are the same
				sigma.append(0)
			else:
				sigma.append(1)                      # the alleles are differen (a mutation and a WT)

		all_sigmas.append(sigma)

tmp_columns = ['name'] + data_positions
sigma_df = pd.DataFrame(all_sigmas, columns=tmp_columns)
sample_names = list(sigma_df.loc[:, 'name'])

# Calculate ds by averaging the correlation matrix
pos_averages = []
for i in range(1, len(tmp_columns)):
	print(i)
	locus_vector = np.array(sigma_df.iloc[:, i])
	numerator = np.nansum(locus_vector)
	denominator = sum(~np.isnan(locus_vector))
	pos_averages.append(float(numerator)/float(denominator))

pos_averages = np.array(pos_averages)
ds = sum(pos_averages)/float(len(pos_averages))
ds_var = np.var(pos_averages)

mcorr_out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/data/mcorr_PSB_data_2024_pipeline_ql_ds_synonymous.csv'
mcorr_out = open(mcorr_out_file, 'w')
mcorr_out.write('# l: the distance between two genomic positions\n')
mcorr_out.write('# m: the mean value of correlatio profile\n')
mcorr_out.write('# v: the variance of correlation profile\n')
mcorr_out.write('# n: the total number of alignments used for calculation\n')
mcorr_out.write('# t: the type of result: Ks is for d_sample, and P2 is for correlation profile\n')
mcorr_out.write('# b: the bootstrap number (all means used all alignments).\n')
mcorr_out.write('l,m,v,n,t,b\n')
mcorr_out.write('0,' + "{:.5f}".format(ds) + ',' + "{:.5f}".format(ds_var) + ',' + str(len(sample_names)) +  ',Ks,all\n')


probabilities_per_distance = {}
for i in range(0, len(data_positions)):
	print(i)
	for j in range(i+1, len(data_positions)):
		pos1, pos2 = data_positions[i], data_positions[j]

		if abs(int(pos1) - int(pos2)) < 600:
			dist = abs(int(pos1) - int(pos2))
			pos1_v, pos2_v = np.array(sigma_df.loc[:, pos1]), np.array(sigma_df.loc[:, pos2])

			# Calculate probability of pos2 having 1 given pos1 has 1  
			# P(B|A) = (P and B) / P(A)
			sum_11_positions = 0.0
			sum_1x_positions = 0.0
			total = 0.0

			for k in range(0, len(pos1_v)):
				print(pos1_v[k], pos2_v[k])
				if pos1_v[k]  in [0, 1] and pos2_v[k] in [0, 1]:
					total += 1
					if pos1_v[k] == 1 and pos2_v[k] == 1:
						sum_11_positions += 1

					if pos1_v[k] == 1:
						sum_1x_positions += 1
					

			if total != 0:
				p_11 = sum_11_positions/total

				if dist not in probabilities_per_distance.keys():
					probabilities_per_distance[dist] = [p_11]
				else:
					probabilities_per_distance[dist].append(p_11)


out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/data/mcorr_PSB_data_2024_ql_ds_synonymous.csv'
out = open(out_file, 'w')

data_out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/data/mcorr_PSB_raw_data_2024_ql_ds_synonymous.csv'
raw_out = open(data_out_file, 'w')
raw_out.write('ds,' + str(ds) + '\n')
# Average over all probabilities per distance (sort it as well)
distances_sorted = list(probabilities_per_distance.keys())
distances_sorted.sort()

distances, probabilities = [], []
for k in distances_sorted:
	distances.append(k)
	v = probabilities_per_distance[k]
	raw_out.write(str(k) + ',' + str(v).strip('[').strip(']') + '\n')
	q_l = sum(v)/len(v)
	probabilities.append(q_l / ds)
	out.write(str(k) + ',' + "{:.5f}".format(sum(v)/len(v)) + '\n')
	np_v = np.array(v)
	var = np.var(np_v)
	mcorr_out.write(str(k) + ',' + "{:.5f}".format(sum(v)/len(v)) + ',' + "{:.5f}".format(var) + ',' + str(len(v)) +',P2,all\n')


out.close()
mcorr_out.close()
raw_out.close()

plt.plot(distances, probabilities, 'o', alpha=0.8, mfc='none')
plt.title('P(l) vs distance for all synonymous PSB loci')
plt.ylabel('average P(l)')
plt.xlabel('Distance l (bp)')
plt.show()
