# four allele test to test how well samples behave as a tree.

from collections import Counter
import numpy as np
import math
import matplotlib.pyplot as plt 
import pandas as pd
import random


def find_pairs(positions, distance, label):
	"""
	Find pairs of loci that are less than a d distance apart from each other.
	We will use them to calculate 4 allele test values
	"""
	good_pairs = []
	counter = 0
	while counter < 300000:
		pos1_i = random.randrange(0, len(positions))
		pos2_i = random.randrange(0, len(positions))
		if abs(int(positions[pos2_i]) - int(positions[pos1_i])) <= distance:
			if pos1_i != pos2_i and (pos1_i, pos2_i) not in good_pairs and (pos2_i, pos1_i) not in good_pairs:
				print(label + ': ' + str(counter))
				good_pairs.append((positions[pos1_i], positions[pos2_i]))
				counter += 1

	return(good_pairs)


def get_allele_freq(positions, good_pairs):
	# This function performs the 4 allele test.
	# It determines two random loci and collects all alleles for the group
	# at that position if neither position is an 'N'.
	# Next, if threshold was satisfied, it checks if there's <= 3 alleles (PASS) or 4 alleles (FAIL)
	# Test looks fine. This function works how it should.

	pass_distance, fail_distance = [], []
	for pair in good_pairs:
		if pair[0] in positions and pair[1] in positions:
		 
			pos1 = list(data.loc[:, pair[0]])
			pos2 = list(data.loc[:, pair[1]])

			list_zip = list(zip(pos1, pos2))
			clean_list = []
			for i in list_zip:
				if not math.isnan(i[0]) and not math.isnan(i[1]):
					clean_list.append(i)

			distance = abs(int(pair[0]) - int(pair[1]))
			if len(clean_list) > 10:
				if len(set(clean_list)) <= 3 and len(set(clean_list)) != 0 :
					#print('Test passed')
					pass_distance.append(distance)
				elif len(set(clean_list)) > 3 :
					#print('Test failed')
					fail_distance.append(distance)
			
	print(len(pass_distance))
	print(len(fail_distance))
	return(pass_distance, fail_distance)


def make_plot_data(pass_d, fail_d, step):
	"""
	For each distance, calculate how many pairs with that distance passed the 4 allele test
	"""
	distances = list(range(0, 50000, step))
	counted_pass = Counter(pass_d)
	counted_fail = Counter(fail_d)

	good_distances = []
	fractions = []
	for i in distances:
		tmp_pass_count, tmp_fail_count = 0, 0
		for j in range(i,i+step):
			if j in counted_pass.keys():
				tmp_pass_count += counted_pass[j]
			if j in counted_fail.keys():
				tmp_fail_count += counted_fail[j]
		if float(tmp_pass_count+tmp_fail_count) != 0:
			fractions.append(tmp_pass_count/float(tmp_pass_count+tmp_fail_count))
			good_distances.append(i)
		

	return(fractions, good_distances)


# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
data_positions = data.columns.tolist()
step = 100

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


good_pairs_full = find_pairs(data_positions, 50000, 'full')
pass_distance, fail_distance = get_allele_freq(data_positions, good_pairs_full)
fractions_full, distances_full = make_plot_data(pass_distance, fail_distance, step)
plt.plot(distances_full, fractions_full, alpha=0.8, label='full data')


good_pairs = find_pairs(positions_5, 50000, '5')
pass_distance, fail_distance = get_allele_freq(data_positions, good_pairs)
fractions_5, distances_5 = make_plot_data(pass_distance, fail_distance, step)
plt.plot(distances_5, fractions_5, alpha=0.8, label='frequency <= 5 SNPs / 1kb')

good_pairs = find_pairs(positions_10, 50000, '10')
pass_distance, fail_distance = get_allele_freq(data_positions, good_pairs)
fractions_10, distances_10 = make_plot_data(pass_distance, fail_distance, step)
plt.plot(distances_10, fractions_10, alpha=0.8, label='frequency <= 10 SNPs / 1kb')

good_pairs = find_pairs(positions_10, 50000, '15')
pass_distance, fail_distance = get_allele_freq(data_positions, good_pairs)
fractions_15, distances_15 = make_plot_data(pass_distance, fail_distance, step)
plt.plot(distances_15, fractions_15, alpha=0.8, label='frequency <= 15 SNPs / 1kb')

plt.ylabel('Fraction of passed 4 allele tests')
plt.xlabel('Distance between loci (bin size = ' + str(step) + 'bp)')
plt.title('Is there a correlation between distance between loci and passing 4 allele test? \n number of pairs=300,000')
plt.legend()
plt.show()

















