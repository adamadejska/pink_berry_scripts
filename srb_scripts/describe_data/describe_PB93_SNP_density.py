##################################################################
# PB93  is a hypermutator. Visualize how many SNPs does it have on its own. 
##################################################################


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
