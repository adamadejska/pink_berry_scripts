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
			pos1 = list(data.iloc[:, pos1_i])[1:]
			pos2 = list(data.iloc[:, pos2_i])[1:]
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


# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_6_no_PB93.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()[1:]
data_positions = data.columns.tolist()

# Count the number of SNPs across all strains in each 1 kb window.
window = 1000
densities = {}
windows = []
for p in data_positions:
    p_int = p.split('.')[0]
    int_div = int(p_int) // window
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
low_density_regions_20 = []
for window, density in densities.items():
    if density <= 10:
        low_density_regions_10.append(window)
    if density <= 5:
        low_density_regions_5.append(window)
    if density <= 15:
        low_density_regions_15.append(window)
    if density <= 20:
        low_density_regions_20.append(window)

# Find the loci that are in the low density regions. We will use them for distance matrix construction.
positions_10, positions_5, positions_15, positions_20 = [], [], [], []
for p in data_positions:
    p_int = p.split('.')[0]
    int_div = int(p_int) // window
    if int_div in low_density_regions_10:
        positions_10.append(p)
    if int_div in low_density_regions_5:
        positions_5.append(p)
    if int_div in low_density_regions_15:
        positions_15.append(p)
    if int_div in low_density_regions_20:
        positions_20.append(p)
        

def find_common_WT_positions():
      #wt_hap_file = "/home/ada/Desktop/PinkBerry_scripts_paper/srb_scripts/haplotyping/haplotype_full_data.txt"
	wt_hap_file = "/home/ada/Desktop/PinkBerry_scripts_paper/srb_scripts/haplotyping/haplotype_full_data_updated_SRB_tree.txt"
	# Create a dictionary where we will store info per contig for each bacteria
	contigs = ['0.0', '5.0', '24.0', '28.0', '67.0', '69.0', '70.0', '87.0']
	contigs_regions = {}

	for i in contigs:
		contigs_regions[i] = []

	# Extract and store the data divided by the strain
	with open(wt_hap_file, 'r') as f:
		for line in f:
			if not line.startswith('['):  # line with contig info and sample
				contig = line.strip().split()[1]
			if line.startswith('['):    # we found a line with wt regions listed
				list_of_regions = list(ast.literal_eval(line))
				wt_positions = []
				for i in list_of_regions:
					wt_positions.extend(list(range(int(i[0]), int(i[1]))))   # Make a list of all possible loci within the WT regions
				
				contigs_regions[contig].append(wt_positions)

	print('figuring out all possible common loci')

	# Find common loci between all strains per contig
	frequent_pos_per_contig = {}
	for c in contigs:
		all_positions_togther = []

		for i in contigs_regions[c]:
			all_positions_togther.extend(i)

		count_positions = Counter(all_positions_togther)
		#print(count_positions)

		frequent_positions = []
		for pos, freq in count_positions.items():
			if freq >= 17:
				frequent_positions.append(pos)
		
		frequent_pos_per_contig[c] = np.array(frequent_positions)

	print('found all possible common loci')

	# Read in the SNP data
	#data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)
	data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_6_no_PB93.csv', index_col=0)

	data_index = data.index.values.tolist()[1:]  # names of bacterial samples
	chr_list = data.iloc[0, :].tolist()
	data_positions = data.columns.tolist()
	data_positions_ints = np.array([int(i.split('.')[0]) for i in data_positions])

	#print(chr_list[:100])

	# Find actual common loci that we have in the biallelic SNPs dataset
	actual_common_loci_per_contig = {}
	for c in contigs:
		contig_mask = [True if i == float(c) else False for i in chr_list]
		contig_positions = data_positions_ints[contig_mask]
		#print(len(contig_positions))
		#print(len(data_positions))
		common_actual_loci = np.intersect1d(contig_positions, frequent_pos_per_contig[c])
		print(len(common_actual_loci))
		actual_common_loci_per_contig[c] = common_actual_loci


	# Only consider positions that are WT and actually in the dataset.
	positions = []
	for i in range(0, len(data_positions)):
		pos = int(data_positions[i].split('.')[0])
		chr = str(chr_list[i])
		if chr in contigs:
			if pos in actual_common_loci_per_contig[chr]:
				positions.append(data_positions[i])

	return(positions)


test_vals_base_10 = get_allele_freq(positions_10, data_index)
test_vals_base_5 = get_allele_freq(positions_5, data_index)
test_vals_base_15 = get_allele_freq(positions_15, data_index)
test_vals_base_20 = get_allele_freq(positions_20, data_index)
test_vals_og = get_allele_freq(data_positions, data_index)
test_vals_base_hap = get_allele_freq(find_common_WT_positions(), data_index)

print('Finished compiling distances')
print('Start making histogram data')
allele_vals_dict_og = Counter(test_vals_og)
allele_vals_dict_base_10 = Counter(test_vals_base_10)
allele_vals_dict_base_5 = Counter(test_vals_base_5)
allele_vals_dict_base_15 = Counter(test_vals_base_15)
allele_vals_dict_base_20 = Counter(test_vals_base_20)
allele_vals_dict_base_hap = Counter(test_vals_base_hap)

N = 4
ind = np.arange(N)  
width = 0.15

def order_values(dict):
	ordered_vals = []
	for i in range(1, 5):
		ordered_vals.append(dict[i])
	return(ordered_vals)


plt.bar(ind, order_values(allele_vals_dict_og), width, label='full data')
plt.bar(ind+width*1, order_values(allele_vals_dict_base_5), width, label='frequency <= 5 SNPs / 1kb')
plt.bar(ind+width*2, order_values(allele_vals_dict_base_10), width, label='frequency <= 10 SNPs / 1kb')
plt.bar(ind+width*3, order_values(allele_vals_dict_base_15), width, label='frequency <= 15 SNPs / 1kb')
plt.bar(ind+width*4, order_values(allele_vals_dict_base_hap), width, label='WT regions based on haplotyping')
#plt.bar(ind+width*3, order_values(allele_vals_dict_base_20), width, label='frequency <= 20 SNPs / 1kb')

plt.title('4 allele test values for SRB loci (coverage >= 6) (n=50,000)')
plt.xlabel('4 allele test values')
plt.ylabel('number of pairs')
plt.yscale('log')
plt.xticks(ind+width,[1, 2, 3, 4])
plt.legend()
plt.show()
