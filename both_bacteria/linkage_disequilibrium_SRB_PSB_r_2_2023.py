################################################################################
# Given a list of positions of two loci and the r^2 value for them,
# create a histogram of all values at certain distances. 
# Linkage disequilibrium analysis 
# Both bacteria PSB and SRB on the same figure.
################################################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# File with all the r^2 data for SRB.
r2_path_srb = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_SRB_data_noPB93_2023_2.csv'

# Initiate bins, each bin will be 10 bp long
bins = list(range(10, 1010, 10))
contigs = [0, 67, 69, 70, 87]
for contig in contigs:
    r_values = {}
    with open(r2_path_srb, 'r') as f:
        for line in f:   # Each line:  locus1, locus2, r^2 value
            locus1, locus2, current_contig, r_2 = line.split(',')
            locus1, locus2 = locus1.split('.')[0], locus2.split('.')[0]
            if int(current_contig) == contig:
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
    plt.plot(bins, mean_r, label='contig ' + str(contig))

# File with all the r^2 data for PSB.
r2_path_psb = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_PSB_data_2023.csv'

# Initiate bins, each bin will be 10 bp long
bins = list(range(0, 1010, 10))

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
mean_r = []
for i in bins:
    mean_r.append(np.mean(np.array(r_values[i])))

# Plot the results
plt.plot(bins, mean_r, label='PSB')
plt.xlabel('Nucleotide Distance (bp) (bin size = 10)')
plt.ylabel('Mean Linkage Disequilibrium r^2 ')
plt.legend()
plt.show()



#####################
# No bins, each distance by itself (bin size = 1)
# All big contigs at once.

# Initiate bins, each bin will be 10 bp long
bins_srb = list(range(5, 1001, 1))
contigs = [0, 67, 69, 70, 87]

for contig in contigs:
    r_values = {}
    with open(r2_path_srb, 'r') as f:
        for line in f:   # Each line:  locus1, locus2, r^2 value
            locus1, locus2, current_contig, r_2 = line.split(',')
            locus1, locus2 = locus1.split('.')[0], locus2.split('.')[0]
            if int(current_contig) == contig:
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
    for i in bins_srb:
        try:
            mean_r.append(np.mean(np.array(r_values[i])))
        except KeyError:
            mean_r.append(0)

    # Plot the results
    plt.plot(bins_srb, mean_r, label='contig ' + str(contig))


### PSB
# Initiate bins, each bin will be 10 bp long
bins_psb = list(range(1, 1001, 1))

r_values = {}
with open(r2_path_psb, 'r') as f:
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
for i in bins_psb:
    try:
        mean_r.append(np.mean(np.array(r_values[i])))
    except KeyError:
        mean_r.append(0)

plt.plot(bins_psb, mean_r, label='PSB')
plt.xlabel('Nucleotide Distance (bp)(log)')
plt.ylabel('Mean Linkage Disequilibrium r^2 (log)')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()

