###############################################################
# Ada Madejska, 2023
# Check SNP frequency along the chromosome of PSB and SRB genomes.
# Are the SNPs clustered into hotspots or is the distribution 
# more constant across the genome? What kind of genes can we find
# in the hotspots? 
###############################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
# Read in the PSB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Count the number of SNPs across all strains in each 1 kb window.
window = 500
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


fig, axs = plt.subplots(2)
fig.suptitle('SNP densities across the genome (PSB)')
axs[0].plot(windows[:1000], list(densities.values())[:1000])
axs[0].set(ylabel='SNPs per kb')
axs[1].plot(windows[1000:2000], list(densities.values())[1000:2000])
axs[1].set(ylabel='SNPs per kb')
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
bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
plt.scatter(bins_mean, n)
plt.title('Histogram of SNPs per 500 bp block (PSB)')
plt.xlabel('Number of SNPs in 500 bp window')
plt.ylabel('count (log scale)')
plt.yscale('log')
plt.show()
"""

### 
# Do similar analysis for SRB but take care of the contig info.
###

# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)
data = data.drop(index=('PB93'))    # PB93 is a hypermutator so get rid of it for this analysis

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
contigs_list = [str(i) for i in set(data.iloc[0,:].tolist())]
#print(data_index)

# [0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
# 20.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 36.0, 37.0, 38.0, 40.0,
# 41.0, 42.0, 43.0, 44.0, 45.0, 49.0, 60.0, 65.0, 66.0, 67.0, 69.0, 70.0, 72.0, 76.0, 77.0, 78.0, 79.0,
# 80.0, 82.0, 83.0, 85.0, 86.0, 87.0, 88.0, 89.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0]

# Count the number of SNPs across all strains in each 1 kb window.
# For now just look at the few biggest contigs: 0, 67, 69, 70, 87
#contigs_list = ['0.0', '67.0', '69.0', '70.0', '87.0']
all_densities = []
all_windows = []
window = 1000
#contigs_list = ['0.0']
#mask = [True if x == 'PB93' else False for x in data_index]
for i in contigs_list:
    densities = {}
    windows = []
    for p in range(0,len(data_positions)):
        if str(data.iloc[0,p]) == i:
            #print(data_positions[p])
            # Check if the position has a mutation in PB93 (supermutant) and ignore any singletons from it.
            #locus_PB93 = np.array(data[data_positions[p]].tolist())[mask]
            locus_all = np.array(data[data_positions[p]].tolist())
            if  sum(locus_all==1) == 0:   # check if there's a mutation in PB93 and it's a singleton 
                continue
            else:
                int_div = int(data_positions[p].split('.')[0]) // window
                if int_div not in densities.keys():
                    densities[int_div] = 1
                    windows.append(int_div)
                else:
                    densities[int_div] += 1

    all_densities.append(densities)
    windows = np.array(windows) / window
    all_windows.append(windows)

# Plot the densities across the genome
fig, axs = plt.subplots(2)
fig.suptitle('SNP densities across the genome (SRB)')
axs[0].plot(all_windows[0][500:800], list(all_densities[0].values())[500:800])
axs[0].set(ylabel='SNPs per kb')
axs[0].title.set_text('Contig 0')
#axs[1].plot(all_windows[contigs_list.index('67.0')][:1000], list(all_densities[contigs_list.index('67.0')].values())[:1000])
#axs[1].set(ylabel='SNPs per kb')
#axs[1].title.set_text('Contig 67')
plt.xlabel('Location in the genome (Mb)')
plt.show()

"""
######################################################################
# Make a histogram for the number of SNPs per kilobase for the whole
# PSB dataset.
######################################################################

# Combine all densities objects
complete_densities = []
for i in all_densities:
    complete_densities = complete_densities + list(i.values())

n, bins, patches = plt.hist(complete_densities, bins=range(1,max(complete_densities)+1, 1))
# Scatter plot
# Find the center of each bin from the bin edges
plt.clf()
bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
plt.scatter(bins_mean, n)
plt.title('Histogram of SNPs per 100 bp block (SRB)')
plt.xlabel('Number of SNPs in 100 bp window')
plt.ylabel('count (log scale)')
plt.yscale('log')
plt.show()

"""
"""
### 
# Do similar analysis for PB93 ONLY since it is a hypermutator.
###

# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
contigs_list = [str(i) for i in set(data.iloc[0,:].tolist())]
#print(data_index)

# Count the number of SNPs across all strains in each 1 kb window.
# For now just look at the few biggest contigs: 0, 67, 69, 70, 87
all_densities = []
all_windows = []
window = 1000
mask = [True if x == 'PB93' else False for x in data_index]
for i in contigs_list:
    densities = {}
    windows = []
    for p in range(0,len(data_positions)):
        if str(data.iloc[0,p]) == i:
            #print(data_positions[p])
            # Check if the position has a mutation in PB93 (supermutant) and ignore any singletons from it.
            locus_PB93 = np.array(data[data_positions[p]].tolist())[mask]
            locus_all = np.array(data[data_positions[p]].tolist())
            if  sum(locus_PB93==1) == 1:   # check if there's a mutation in PB93  
                int_div = int(data_positions[p].split('.')[0]) // window
                if int_div not in densities.keys():
                    densities[int_div] = 1
                    windows.append(int_div)
                else:
                    densities[int_div] += 1

                
    all_densities.append(densities)
    windows = np.array(windows) / window
    all_windows.append(windows)

# Plot the densities across the genome
fig, axs = plt.subplots(2)
fig.suptitle('SNP densities across the genome (SRB) PB93 ONLY')
axs[0].plot(all_windows[0][:1000], list(all_densities[0].values())[:1000])
axs[0].set(ylabel='SNPs per kb')
axs[0].title.set_text('Contig 0')
axs[1].plot(all_windows[contigs_list.index('67.0')][:1000], list(all_densities[contigs_list.index('67.0')].values())[:1000])
axs[1].set(ylabel='SNPs per kb')
axs[1].title.set_text('Contig 67')
plt.xlabel('Location in the genome (Mb)')
plt.show()
"""