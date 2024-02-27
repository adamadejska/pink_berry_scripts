################################################################################
# Analysis similar to Sakoparing et al. "Whole genome phylogenies reflect the distribution of recombination rates [...]"
# Given a pair of strains, find which SNP loci are different between those two strains.
# i.e If strain_1 has a mutation at pos1 and strain_2 does not, or vice versa, count that as a mutation.
# Create a SNP density graph based on those mutations.
# See if there are any regions with high SNP density. (possible HGT events).
# Also note the overall pairwise divergence for each pair - we wanna see how SNP density graph
# changes based on the pairwise divergence.
################################################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# PSB, read the data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
"""
diversities = {}
for strain1_i in range(0, len(data_index)):
    for strain2_i in range(strain1_i+1, len(data_index)):
        
        strain1 = np.array(data.iloc[strain1_i, :])
        strain2 = np.array(data.iloc[strain2_i, :])

        # diversity = number_of_differences / length of strain
        diversity_count = 0
        length = 0
        for i in range(0, len(strain1)):
            if not np.isnan(strain1[i]) and not np.isnan(strain2[i]):
                length += 1
                if strain1[i] != strain2[i]:
                    diversity_count += 1

        diversity = diversity_count / length
        diversities[(data_index[strain1_i], data_index[strain2_i])] = diversity
"""
# try pairs:
# ('PB28', 'PB73'): 0.21915908812854287  -> high diversity
# ('AAGAGGCA-ACTGCATA', 'PB16'): 0.07362483565884897    -> medium diversity
# ('PB24', 'PB65'): 0.07655255983035444  -> medium diversity
# ('PB_1', 'PB25'): 0.0033328976604365445     -> low diversity
        
#print({k: v for k, v in sorted(diversities.items(), key=lambda item: item[1])})
        
strain1_name = 'PB_1'
strain2_name = 'PB25'

strain1 = np.array(data.loc[strain1_name, :])
strain2 = np.array(data.loc[strain2_name, :])

diverse_positions = []
for i in range(0, len(data_positions)):
    if not np.isnan(strain1[i]) and not np.isnan(strain2[i]):
        if strain1[i] != strain2[i]:
            diverse_positions.append(data_positions[i])

print(len(diverse_positions))
# Count the number of SNPs across all strains in each 1 kb window.
window = 1000
densities = {}
windows = list(range(0, 8000))
for i in windows:
    densities[i] = 0

for p in diverse_positions:
    int_div = int(p) // window
    densities[int_div] += 1

# Plot the densities across the genome
windows = np.array(windows) / window


plt.title('SNP densities across the genome (PSB)')
plt.plot(windows[:2500], list(densities.values())[:2500])
plt.ylabel('SNPs per kb')
plt.xlabel('Location in the genome (Mb)')
plt.show()


######################################################################
# Make a histogram for the number of SNPs per kilobase for the whole
# PSB dataset.
######################################################################

n, bins, patches = plt.hist(list(densities.values()), bins=range(1,max(densities.values())+1, 1), edgecolor='black')
# Scatter plot
# Find the center of each bin from the bin edges
plt.clf()
#bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
#plt.scatter(bins_mean, n)
plt.scatter(bins[:-1], n)
plt.title('Histogram of SNPs per 1000 bp block (PSB)')
plt.xlabel('Number of SNPs in 1000 bp window')
plt.ylabel('count (log scale)')
plt.yscale('log')
plt.show()