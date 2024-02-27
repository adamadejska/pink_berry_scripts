################################################################################
# Given a list of positions of two loci and the r^2 value for them,
# create a histogram of all values at certain distances. 
# Linkage disequilibrium analysis 
# Both bacteria PSB and SRB on the same figure.
################################################################################


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import numpy as np
import pandas as pd

# File with all the r^2 data for SRB.
r2_path_srb = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_SRB_data_noPB93_2023_2.csv'

# Initiate bins, each bin will be 10 bp long
bins_srb = list(range(10, 1010, 10))
#contigs = [0, 67, 69, 70, 87]
#for contig in contigs:
r_values = {}
with open(r2_path_srb, 'r') as f:
    for line in f:   # Each line:  locus1, locus2, r^2 value
        locus1, locus2, current_contig, r_2 = line.split(',')
        locus1, locus2 = locus1.split('.')[0], locus2.split('.')[0]
            #if int(current_contig) == contig:
        if locus1 == locus2:
            continue
        else:
            dist = int(locus2) - int(locus1)
            rounded_down = dist - dist % 10
            
            if rounded_down not in r_values.keys():
                r_values[rounded_down] = [float(r_2)]
            else:
                r_values[rounded_down].append(float(r_2))

# Calculate average values for each bin.
mean_r_srb = []
for i in bins_srb:
    mean_r_srb.append(np.mean(np.array(r_values[i])))


# File with all the r^2 data for PSB.
r2_path_psb = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_PSB_data_2023.csv'

# Initiate bins, each bin will be 10 bp long
bins_psb = list(range(10, 1010, 10))

r_values = {}
with open(r2_path_psb, 'r') as f:
    for line in f:   # Each line:  locus1, locus2, r^2 value
        locus1, locus2, r_2 = line.split(',')
        if locus1 == locus2:
            continue
        else:
            dist = int(locus2) - int(locus1)
            rounded_down = dist - dist % 10
            
            if rounded_down not in r_values.keys():
                r_values[rounded_down] = [float(r_2)]
            else:
                r_values[rounded_down].append(float(r_2))

# Calculate average values for each bin.
mean_r_psb = []
for i in bins_psb:
    mean_r_psb.append(np.mean(np.array(r_values[i])))


# Find a best fit line for both psb and srb
np_bins_srb = np.array(bins_srb)
np_mean_r_srb = np.array(mean_r_srb)
log_mean_r_srb = np.log(np.array(mean_r_srb))
a_srb, b_srb = np.polyfit(np_bins_srb, np_mean_r_srb, 1)
y = 0.9 / pow(np_bins_srb, 0.3)
print(a_srb)
print(b_srb)

# Plot the results
fig, ax1 = plt.subplots()
ax1.plot(bins_srb, mean_r_srb, label='SRB', color='blue')
ax1.plot(bins_psb, mean_r_psb, label='PSB', color='orange')
ax1.plot(np_bins_srb, y, c='r')
plt.title('PSB and SRB $r^2$ measurements')
ax1.set_xlabel('Nucleotide Distance (bp) (bin size = 10) (log)')
ax1.set_ylabel('Mean Linkage Disequilibrium $r^2$ (log)')
ax1.legend()
ax1.set_yscale('log')
ax1.set_xscale('log')
plt.show()
"""

# this is an inset axes over the main axes
ax2 = plt.axes([1,1,0,0])
ip = InsetPosition(ax1, [0.1,0.1,0.35,0.35])
ax2.set_axes_locator(ip)

r2_path = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_PSB_data_2023.csv'

# Initiate bins, each bin will be 10 bp long
bins = list(range(10, 1001, 10))

r_values = {}
with open(r2_path, 'r') as f:
    for line in f:   # Each line:  locus1, locus2, r^2 value
        locus1, locus2, r_2 = line.split(',')
        if locus1 == locus2:
            continue
        else:
            dist = int(locus2) - int(locus1)
            
            if dist not in r_values.keys():
                r_values[dist] = [float(r_2)]
            else:
                r_values[dist].append(float(r_2))

# Calculate average values for each bin.
mean_r = []
for i in bins:
    try:
        mean_r.append(np.mean(np.array(r_values[i])))
    except KeyError:
        mean_r.append(0)

ax2.plot(bins, mean_r, color='orange')

# File with all the r^2 data.
r2_path = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_SRB_data_noPB93_2023_2.csv'

# Initiate bins, each bin will be 10 bp long
bins = list(range(10, 1010, 10))
r_values = {}
with open(r2_path, 'r') as f:
    for line in f:   # Each line:  locus1, locus2, r^2 value
        locus1, locus2, current_contig, r_2 = line.split(',')
        locus1, locus2 = locus1.split('.')[0], locus2.split('.')[0]
        if locus1 == locus2:
            continue
        else:
            dist = int(locus2) - int(locus1)
            rounded_down = dist - dist % 10
            
            if rounded_down not in r_values.keys():
                r_values[rounded_down] = [float(r_2)]
            else:
                r_values[rounded_down].append(float(r_2))

# Calculate average values for each bin.
mean_r = []
for i in bins:
    mean_r.append(np.mean(np.array(r_values[i])))

# Plot the results
ax2.plot(bins, mean_r, color='blue')

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel('Nucleotide Distance (bp)')
ax2.set_ylabel('Mean Linkage Disequilibrium $r^2$')
#ax2.set_title('PSB r^2')
plt.show()
"""