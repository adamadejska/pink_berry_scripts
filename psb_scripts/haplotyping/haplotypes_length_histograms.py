"""
Make a histogram of lengths for green haplotype (PSB). 
"""
import matplotlib.pyplot as plt
import math
import numpy as np


lengths_file_path = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/haplotypes_lengths.txt'


data_dict = {}
# Read in the lengths file for different conditions
with open(lengths_file_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            condition = line[1:].strip()
            green_lengths = f.readline().replace('[','').replace(']','') # get rid of any unnecessary characters
            green_lengths = [int(i) for i in green_lengths.split(',')] # change numbers to ints
            green_lengths = [i for i in green_lengths if i != 0]   # clean zeros

            data_dict[condition] = green_lengths


######################################################################################
# Show a cumulant graph for the distribution of lengths given different conditions

bins = np.arange(0, 125000, 1000)
cumulative_data = {}
for condition in data_dict.keys():
    sizes, bins, _ = plt.hist(data_dict[condition], bins=bins)

    total = float(sum(sizes))

    # Calculate the cumulant for minimap2 identities
    fractions = [0.0]
    current_total = 0.0
    for i in sizes:
        current_total += (i/total)
        fractions.append(current_total)

    cumulative_data[condition] = fractions

plt.clf()

for condition in cumulative_data.keys():
    plt.plot(bins, cumulative_data[condition], label=condition)


"""
for condition in data_dict.keys():
    # Make a histogram for each haplotype separately
    bins = range(0, max(data_dict[condition][0])+1000, 1000)
    plt.hist(data_dict[condition][0], bins=bins, label=condition)
"""
plt.xlabel('haplotype block length (bp)  (bin size=1000bp)')
plt.ylabel('Fraction of blocks')
plt.title('Cummulative graph of the green haplotype blocks frequency by length \n all PSB strains')
plt.legend()
plt.show()

"""
# Plot a histogram of just one condition we focus on: forward, max penalty=5, L=5
lengths_green = data_dict['forward (max penalty = 5, L=5)']
bins = np.arange(0, max(lengths_green)+1000, 1000)



plt.clf()
y, x, _ = plt.hist(lengths_green, bins=bins, edgecolor='black')
x_full = x.copy()
y = np.array([math.log(i) for i in y[:30]])
x = np.array(x[:30])
# Best fit line
a, b = np.polyfit(x, y, 1)

plt.clf()
  
plt.hist(lengths_green, bins=bins, edgecolor='black')
plt.title('Green haplotype blocks frequency by length\n mixing layer strains only, max penalty=5, L=5')
plt.xlabel('haplotype block length (bp) (bin size=1000bp)')
plt.ylabel('Number of blocks (log)')
plt.yscale('log')
x = x_full[:55]
plt.plot(x, np.exp(a*x+b))
plt.show()
"""