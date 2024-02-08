# four allele test to test how well samples behave as a tree.

from collections import Counter
import numpy as np
import math
import matplotlib.pyplot as plt 
import pandas as pd
import random

# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_6_no_PB93.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()[1:]
positions = data.columns.tolist()[1:]


def get_allele_freq(positions, index):
	# This function performs the 4 allele test.
	# It determines two random loci and collects all alleles for the group
	# at that position if neither position is an 'N'.
	# Next, if threshold was satisfied, it checks if there's <= 3 alleles (PASS) or 4 alleles (FAIL)

	test_vals = []
	counter = 0
	pairs_used = []
	while counter <= 500000:
		print(counter)
		# Get random positions
		pos1_i = random.randrange(0, len(positions))
		pos2_i = random.randrange(0, len(positions))

		if pos1_i != pos2_i and (pos1_i, pos2_i) not in pairs_used and (pos2_i, pos1_i) not in pairs_used:
			pos1 = list(data.iloc[:, pos1_i])[1:]
			pos2 = list(data.iloc[:, pos2_i])[1:]

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


def make_hist_data(pass_d, fail_d, max_len):
	# This function makes a pandas dataframe that will be used as a histogram input.
	# It separated the data based on bins.

	# Sort the data based on bins.
	d_len = {}
	for i in range(0, max_len, 50000):
		d_len[i] = [[], []]
	#d_len[bin] = [[pass], [fail]]

	for i in pass_d:
		for key in d_len.keys():
			if i >= key and i < key+50000:
				d_len[key][0].append(i)
	for i in fail_d:
		for key in d_len.keys():
			if i >= key and i < key+50000:
				d_len[key][1].append(i)

	keys = d_len.keys()
	keys.sort()
	keys = [i/1000 for i in keys]

	# Calculate the fractions for the bar graph.
	cdf = pd.DataFrame(0.0, index=keys, columns=['pass'])
	for k, v in d_len.items():
		#print(v)
		key = k / 1000
		total = float(len(v[0]) + len(v[1]))
		if total > 0:
			pass_per = float(len(v[0])) / total
			cdf.at[key, 'pass'] = pass_per * 100
		else:
			cdf.at[key, 'pass'] = 0

	return(cdf)


def plot_hist(cdf, name):
# This function plots the histogram.

	
	print('Finished making histogram data.')
	print('Start plotting data.')
	cdf.plot(kind='bar')
	ticks = cdf.index.values.tolist()
	for i in range(0, len(ticks)):
		if i % 5 != 0:
			ticks[i] = ''
	plt.xticks(range(len(ticks)), ticks, rotation='vertical')
	plt.xlabel('Length between two loci (Kb)')
	plt.ylabel('Percentage of failed tests (%)')
	plt.title('Group ' + name + ' probability of passing 4 allele test based on distance between loci. \n n = 50,000 \n threshold n >= 6' )
	plt.show()

####################################
# MAIN
	
test_vals = get_allele_freq(positions, index)
print('Finished compiling distances')
print('Start making histogram data')
allele_vals_dict = Counter(test_vals)
plt.bar(allele_vals_dict.keys(), allele_vals_dict.values())
plt.title('4 allele test values for SRB loci (full data, n=500,000)')
plt.xlabel('4 allele test values')
plt.ylabel('number of pairs')
plt.yscale('log')
plt.xticks([1, 2, 3, 4])
plt.show()
