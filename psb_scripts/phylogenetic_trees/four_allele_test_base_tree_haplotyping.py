import ast
import numpy as np
import pandas as pd
from collections import Counter
import random
import math
import matplotlib.pyplot as plt

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
	if freq >= num_of_bacteria-5:
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
data_positions = [str(i).split('.')[0] for i in data_positions.tolist()]


def get_allele_freq(positions, index):
	# This function performs the 4 allele test.
	# It determines two random loci and collects all alleles for the group
	# at that position if neither position is an 'N'.
	# Next, if threshold was satisfied, it checks if there's <= 3 alleles (PASS) or 4 alleles (FAIL)

	test_vals = []
	counter = 0
	pairs_used = []
	while counter <= 50000:
		print(counter)
		# Get random positions
		pos1_i = random.randrange(0, len(positions))
		pos2_i = random.randrange(0, len(positions))

		if pos1_i != pos2_i and (pos1_i, pos2_i) not in pairs_used and (pos2_i, pos1_i) not in pairs_used:
			pos1 = list(data.iloc[:, pos1_i])
			pos2 = list(data.iloc[:, pos2_i])
			allele = set()

			list_zip = list(zip(pos1, pos2))
			clean_list = []
			for i in list_zip:
				if not math.isnan(i[0]) and not math.isnan(i[1]):
					clean_list.append(i)
			
			if len(clean_list) > 10:
				test_vals.append(len(set(clean_list)))
				counter += 1
				pairs_used.append((pos1_i, pos2_i))

	return(test_vals)

test_vals_base = get_allele_freq(common_positions, data_index)
allele_vals_dict_base = Counter(test_vals_base)
N = 4
ind = np.arange(N)  
width = 0.2
plt.bar(ind, allele_vals_dict_base.values(), width, label='frequency <= 5 SNPs / 1kb')
plt.show()