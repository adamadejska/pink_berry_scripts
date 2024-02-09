###############################################################
# Ada Madejska, 2023
# Create a base tree for PSB. Given our analysis of SNP frequency 
# across the genome, create a tree based on only the loci that
# were in low frequency regions in that previous analysis.
# We want to establish some kind of base using well conserved regions.  
###############################################################

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
data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_3_no_PB93.csv', index_col=0)

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

test_vals_base_10 = get_allele_freq(positions_10, data_index)
#test_vals_base_5 = get_allele_freq(positions_5, data_index)
test_vals_base_15 = get_allele_freq(positions_15, data_index)
test_vals_base_20 = get_allele_freq(positions_20, data_index)
test_vals_og = get_allele_freq(data_positions, data_index)

print('Finished compiling distances')
print('Start making histogram data')
allele_vals_dict_og = Counter(test_vals_og)
allele_vals_dict_base_10 = Counter(test_vals_base_10)
#allele_vals_dict_base_5 = Counter(test_vals_base_5)
allele_vals_dict_base_15 = Counter(test_vals_base_15)
allele_vals_dict_base_20 = Counter(test_vals_base_20)

N = 4
ind = np.arange(N)  
width = 0.2
plt.bar(ind, allele_vals_dict_og.values(), width, label='full data')
#plt.bar(ind+width, allele_vals_dict_base_5.values(), width, label='frequency <= 5 SNPs / 1kb')
plt.bar(ind+width*1, allele_vals_dict_base_10.values(), width, label='frequency <= 10 SNPs / 1kb')
plt.bar(ind+width*2, allele_vals_dict_base_15.values(), width, label='frequency <= 15 SNPs / 1kb')
plt.bar(ind+width*3, allele_vals_dict_base_20.values(), width, label='frequency <= 20 SNPs / 1kb')

plt.title('4 allele test values for SRB loci (n=50,000)')
plt.xlabel('4 allele test values')
plt.ylabel('number of pairs')
plt.yscale('log')
plt.xticks(ind+width,[1, 2, 3, 4])
plt.legend()
plt.show()
