################################################################################
# Given a list of positions of two loci and the r^2 value for them,
# create a histogram of all values at certain distances. 
# Linkage disequilibrium analysis SRB 2023
################################################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# File with all the r^2 data.
r2_path = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_SRB_data_noPB93_2023_2.csv'

# Initiate bins, each bin will be 10 bp long
bins = list(range(10, 1010, 10))
contigs = [0, 67, 69, 70, 87]
for contig in contigs:
    r_values = {}
    with open(r2_path, 'r') as f:
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

plt.title('SRB r^2 measurements')
plt.xlabel('Nucleotide Distance (bp) (bin size = 10) ')
plt.ylabel('Mean Linkage Disequilibrium r^2 ')
plt.legend()
plt.show()

#####################
# No bins, each distance by itself (bin size = 1)

# Initiate bins, each bin will be 10 bp long
bins = list(range(5, 1001, 1))

contig = 0
r_values = {}
with open(r2_path, 'r') as f:
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
for i in bins:
    try:
        mean_r.append(np.mean(np.array(r_values[i])))
    except KeyError:
        mean_r.append(0)

#print(bins)
# Plot the results
plt.plot(bins, mean_r)
plt.title('SRB r^2 measurements (contig 0)')
plt.xlabel('Nucleotide Distance (bp)')
plt.ylabel('Mean Linkage Disequilibrium r^2')
plt.show()

# Plot the results as log-log 
plt.plot(bins, mean_r)
plt.title('SRB r^2 measurements (contig 0)')
plt.xlabel('Nucleotide Distance (bp) (log scale)')
plt.ylabel('Mean Linkage Disequilibrium r^2 (log scale)')
plt.yscale('log')
plt.xscale('log')
plt.show()

#####################
# No bins, each distance by itself (bin size = 1)
# All big contigs at once.

# Initiate bins, each bin will be 10 bp long
bins = list(range(5, 1001, 1))
contigs = [0, 67, 69, 70, 87]

for contig in contigs:
    r_values = {}
    with open(r2_path, 'r') as f:
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
    for i in bins:
        try:
            mean_r.append(np.mean(np.array(r_values[i])))
        except KeyError:
            mean_r.append(0)

    # Plot the results
    plt.plot(bins, mean_r, label='contig ' + str(contig))

plt.title('SRB r^2 measurements')
plt.xlabel('Nucleotide Distance (bp)(log)')
plt.ylabel('Mean Linkage Disequilibrium r^2 (log)')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()

