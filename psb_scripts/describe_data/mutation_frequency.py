####################################################################
# What's the mutation frequency per locus? 
# How many singletons / doubletons / etc ?
####################################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in the PSB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
#data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/psb/PSB_snp_data_coverage_6.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

mutation_counts = []
for p in data_positions:
    tmp_loc = np.array(data.loc[:, p])
    mutation_counts.append(sum(tmp_loc==1))

# Create a histogram
plt.hist(mutation_counts, bins=range(0, max(mutation_counts)+1, 1), edgecolor='black')
plt.yscale('log')
plt.xlabel('number of mutations per locus')
plt.ylabel('number of loci (log scale)')
plt.title('Mutation frequency PSB')
plt.show()