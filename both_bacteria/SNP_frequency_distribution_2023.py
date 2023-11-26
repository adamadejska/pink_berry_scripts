###############################################################
# Ada Madejska, 2023
# What is the frequency of SNPs in PSB and SRB?
###############################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# PSB, read the data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

frequencies_psb = []
for p in data_positions:
    locus = np.array(data.loc[:, p])
    frequencies_psb.append(sum(locus==1))

################################################################
# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

frequencies_srb = []
for p in data_positions:
    locus = np.array(data.loc[:, p])[1:]
    frequencies_srb.append(sum(locus==1))
"""
# plot frequency distribution
fig, axes = plt.subplots(nrows=1, ncols=2)
axes[0].hist(frequencies_psb, bins=range(0, max(frequencies_psb)+1, 1),edgecolor='black')
axes[1].hist(frequencies_srb, bins=range(0, max(frequencies_srb)+1, 1),edgecolor='black')

axes[0].set_title('PSB')
axes[1].set_title('SRB')

axes[0].set_yscale('log')
axes[1].set_yscale('log')

axes[0].set_xscale('log')
axes[1].set_xscale('log')

axes[0].set_ylabel('frequency (log)')
fig.supxlabel('number of SNPs per locus')
fig.suptitle('SNP frequency specturm')
plt.show()
"""
n1, bins1, patches1 = plt.hist(frequencies_psb, bins=range(0, max(frequencies_psb)+1, 1),edgecolor='black')
n2, bins2, patches2 = plt.hist(frequencies_srb, bins=range(0, max(frequencies_srb)+1, 1),edgecolor='black')
plt.clf()
# Scatter plot
# Find the center of each bin from the bin edges
fig, axes = plt.subplots(nrows=1, ncols=2)
bins_mean1 = [0.5 * (bins1[i] + bins1[i+1]) for i in range(len(n1))]
axes[0].scatter(bins_mean1, n1, label='data')

bins_mean2 = [0.5 * (bins2[i] + bins2[i+1]) for i in range(len(n2))]
axes[1].scatter(bins_mean2, n2, label='data')

# Add known slopes of -1 and -2 to the graphs to compare to the experimental data.
slope1_x = np.array([(1/i ) for i in range(2,30)])   # slope = -1 on log log
print(slope1_x[0] - slope1_x[-1] )
slope2_x = np.array([(1/pow(i,2) * 80000) for i in range(9,100)])   # slope = -2 on log log
slope3_x = np.array([(1/pow(i,3) * 900000) for i in range(10,100)])   # slope = -3 on log log
y_vals1 = np.array(range(2,30))
y_vals2 = np.array(range(9,100))
y_vals3 = np.array(range(10,100))
axes[0].plot(y_vals1, slope1_x, c='r', label='slope= -1')
axes[0].plot(y_vals2, slope2_x, c='g', label='slope= -2')
axes[0].plot(y_vals3, slope3_x, c='m', label='slope= -3')

slope1_x = np.array([(20000/i) for i in range(2,50)])   # slope = -1 on log log
slope2_x = np.array([(200000/pow(i,2)) for i in range(10,80)])   # slope = -2 on log log
slope3_x = np.array([(7000000/pow(i,3)) for i in range(25,100)])   # slope = -3 on log log
y_vals1 = np.array(range(2,50))
y_vals2 = np.array(range(10,80))
y_vals3 = np.array(range(25,100))
axes[1].plot(y_vals1, slope1_x, c='r', label='slope= -1')
axes[1].plot(y_vals2, slope2_x, c='g', label='slope= -2')
axes[1].plot(y_vals3, slope3_x, c='m', label='slope= -3')

axes[0].set_title('PSB')
axes[1].set_title('SRB')

axes[0].set_yscale('log')
axes[1].set_yscale('log')

axes[0].set_xscale('log')
axes[1].set_xscale('log')

axes[0].legend() 
axes[1].legend() 

axes[0].set_ylabel('frequency (log)')
fig.supxlabel('number of SNPs per locus (log)')
#fig.suptitle('SNP frequency specturm')
plt.show()