###############################################################
# Ada Madejska, 2023
# Create a base tree for PSB. Given our analysis of SNP frequency 
# across the genome, create a tree based on only the loci that
# were in low frequency regions in that previous analysis.
# We want to establish some kind of base using well conserved regions.  
###############################################################

import ast
from collections import Counter
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import random


def find_pairs(positions, distance):
	"""
	Find pairs of loci that are less than a d distance apart from each other.
	We will use them to calculate 4 allele test values
	"""
	good_pairs = []
	counter = 0
	while counter < 50000:
		pos1_i = random.randrange(0, len(positions))
		pos2_i = random.randrange(0, len(positions))
		if abs(int(positions[pos2_i]) - int(positions[pos1_i])) >= distance:
			if pos1_i != pos2_i and (pos1_i, pos2_i) not in good_pairs and (pos2_i, pos1_i) not in good_pairs:
				print(counter)
				good_pairs.append((positions[pos1_i], positions[pos2_i]))
				counter += 1

	return(good_pairs)


def get_allele_freq(positions, good_pairs, clade):
	# This function performs the 4 allele test.
	# It determines two random loci and collects all alleles for the group
	# at that position if neither position is an 'N'.
	# Next, if threshold was satisfied, it checks if there's <= 3 alleles (PASS) or 4 alleles (FAIL)

	test_vals = []
	for pair in good_pairs:
		if pair[0] in positions and pair[1] in positions:
		 
			pos1 = list(data.loc[clade, pair[0]])
			pos2 = list(data.loc[clade, pair[1]])

			list_zip = list(zip(pos1, pos2))
			clean_list = []
			for i in list_zip:
				if not math.isnan(i[0]) and not math.isnan(i[1]):
					clean_list.append(i)
			
			
			test_vals.append(len(set(clean_list)))

	return(test_vals)


def find_common_WT_positions(n):
	 # Read in the wt haplotype file
	wt_hap_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/haplotyping/haplotypes_full_data.txt'

	print('start')
	# Create a numpy matrix to store the data divided by the strain
	m = []
	num_of_bacteria = 0
	with open(wt_hap_file, 'r') as f:
		for line in f:
			if not line.startswith('#'):    # we found a line with wt regions listed
				#line = line[1:-1]
				list_of_regions = list(ast.literal_eval(line))
				wt_positions = np.array([])
				for i in list_of_regions:
					wt_positions = np.append(wt_positions, range(int(i[0]), int(i[1])))   # Make a list of all possible loci within the WT regions
				
				m.append(wt_positions)
				num_of_bacteria += 1

	print(num_of_bacteria)
	print('figuring out all possible common loci')
	# Find common loci between all strains

	all_positions_togther = []

	for i in m:
		all_positions_togther.extend(i)

	count_positions = Counter(all_positions_togther)
	#print(count_positions)

	frequent_positions = []
	for pos, freq in count_positions.items():
		if freq >= num_of_bacteria-n:   # how strict do we want to be with choosing positions?
			frequent_positions.append(pos)

	common = np.array(frequent_positions)

	print('found all common possible loci')
	# Read in the SNP data
	data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
	clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']

	data_index = data.index.values.tolist()  # names of bacterial samples
	data_positions = np.array([int(i) for i in data.columns.tolist()])

	common_actual_loci = np.intersect1d(data_positions, common)

	print(len(common_actual_loci))

	# Only consider positions that are not singletons and that are not fixed in our current test clade.
	common_positions = [str(i).split('.')[0] for i in common_actual_loci.tolist()]

	return(common_positions)

# Read in the PSB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
#data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/PSB_snp_data_coverage_6.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Count the number of SNPs across all strains in each 1 kb window.
window = 1000
densities = {}
windows = []
for p in data_positions:
	int_div = int(p) // window
	if int_div not in densities.keys():
		densities[int_div] = 1
		windows.append(int_div)
	else:
		densities[int_div] += 1

# Plot the densities across the genome
windows = np.array(windows) / window

# Filter the windows based on the number of mutations they contain
# We want only the windows that have less than 10 SNPs per 1 kb window (for PSB)

low_density_regions_10 = []
low_density_regions_5 = []
low_density_regions_15 = []
for window, density in densities.items():
	if density <= 10:
		low_density_regions_10.append(window)
	if density <= 5:
		low_density_regions_5.append(window)
	if density <= 15:
		low_density_regions_15.append(window)

# Find the loci that are in the low density regions. We will use them for distance matrix construction.
positions_10, positions_5, positions_15 = [], [], []
for p in data_positions:
	int_div = int(p) // window
	if int_div in low_density_regions_10:
		positions_10.append(p)
	if int_div in low_density_regions_5:
		positions_5.append(p)
	if int_div in low_density_regions_15:
		positions_15.append(p)

clade_9_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']
ladder = ['GGACTCCT-CTAAGCCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB_8','PB24','PB64','PB55','AAGAGGCA-AGAGTAGA','GCTACGCT-GTAAGGAG','AAGAGGCA-ACTGCATA','AAGAGGCA-AAGGAGTA','GCTACGCT-ACTGCATA','PB31','AAGAGGCA-GTAAGGAG','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB32','PB16','PB47','PB48']
bc_strains = ['AGGCAGAA-ACTGCATA','AGGCAGAA-CTAAGCCT','TAAGGCGA-TATCCTCT','PB78','PB61','GGACTCCT-TAGATCGC','PB45','PB77','PB69','PB_5','PB37','PB13','PB85','PB28','PB67','PB59','PB52','PB26','PB34','PB11']
e_clade = ['PB11','CTCTCTAC-AAGGAGTA','GTAGAGGA-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','CTCTCTAC-TAGATCGC','CTCTCTAC-ACTGCATA','CAGAGAGG-TAGATCGC','CAGAGAGG-CTAAGCCT','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','CAGAGAGG-ACTGCATA','CAGAGAGG-AGAGTAGA','CAGAGAGG-TATCCTCT','CAGAGAGG-AAGGAGTA','CTCTCTAC-AGAGTAGA','CTCTCTAC-GTAAGGAG','PB34','PB26']
basal_clade = ['PB29','CGTACTAG-TATCCTCT','PB57','PB51','AGGCAGAA-GTAAGGAG','PB90','PB36','PB82','PB20','PB65','PB60','PB21','PB58','PB43','PB66','GTAGAGGA-AAGGAGTA','PB12','AGGCAGAA-AAGGAGTA','GTAGAGGA-CTCTCTAT','PB75','AGGCAGAA-AGAGTAGA','PB41','CGTACTAG-TAGATCGC','PB25','PB_2']

clades = [clade_9_strains, ladder, bc_strains, e_clade, basal_clade]
clade_name = ['9 clade', 'ladder', 'BC clade', 'E clade', 'basal clade']
good_pairs = find_pairs(data_positions, 0)

for k in range(0, len(clades)):
	clade = clades[k]
	test_vals_og = get_allele_freq(data_positions, good_pairs, clade)
	test_vals_base_10 = get_allele_freq(positions_10, good_pairs, clade)
	test_vals_base_5 = get_allele_freq(positions_5, good_pairs, clade)
	test_vals_base_15 = get_allele_freq(positions_15, good_pairs, clade)
	#test_vals_base_hap_5 = get_allele_freq(find_common_WT_positions(5))
	#test_vals_base_hap_strict = get_allele_freq(find_common_WT_positions(0))


	print('Finished compiling distances')
	print('Start making histogram data')
	allele_vals_dict_og = Counter(test_vals_og)
	#allele_vals_dict_base_hap_5 = Counter(test_vals_base_hap_5)
	#allele_vals_dict_base_hap_strict = Counter(test_vals_base_hap_strict)
	allele_vals_dict_base_10 = Counter(test_vals_base_10)
	allele_vals_dict_base_5 = Counter(test_vals_base_5)
	allele_vals_dict_base_15 = Counter(test_vals_base_15)

	N = 4
	ind = np.arange(N)  
	width = 0.1

	def order_values(dict):

		ordered_vals = []

		# Also normalize the values
		total = 0
		for i in range(1, 5):
			ordered_vals.append(dict[i])
			total += dict[i]

		normalized_vals = list(np.array(ordered_vals)/total)

		return(normalized_vals)

	plt.bar(ind, order_values(allele_vals_dict_og), width, label='full data')
	plt.bar(ind+width, order_values(allele_vals_dict_base_5), width, label='frequency <= 5 SNPs / 1kb')
	plt.bar(ind+width*2, order_values(allele_vals_dict_base_10), width, label='frequency <= 10 SNPs / 1kb')
	plt.bar(ind+width*3, order_values(allele_vals_dict_base_15), width, label='frequency <= 15 SNPs / 1kb')
	#plt.bar(ind+width*4, order_values(allele_vals_dict_base_hap_5), width, label='WT regions based on haplotyping (90%)')
	#plt.bar(ind+width*5, order_values(allele_vals_dict_base_hap_strict), width, label='WT regions based on haplotyping (strict)')

	plt.title('4 allele test values for PSB loci (coverage >= 3) (n=50,000)\n' + clade_name[k])
	plt.xlabel('4 allele test values')
	plt.ylabel('fraction of pairs (normalized)')
	plt.yscale('log')
	plt.xticks(ind+width*2,[1, 2, 3, 4])
	plt.legend()
	plt.show()
