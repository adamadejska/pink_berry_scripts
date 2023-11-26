###########################################################################
# Make a plot with three haplotypes graphs for different variants of the variables
# PSB
###########################################################################

import matplotlib.pyplot as plt

haplotype_file = '/home/ada/Desktop/Shraiman_lab/9_clade_analysis/haplotyping/AAGAGGCA_haplotyping.txt'

conditions = {}
with open(haplotype_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            line = line.split()
            chr, condition = line[2], ' '.join(line[3:])
            haplotype = [int(i) for i in f.readline().strip().split(',')]
            if chr not in conditions.keys():
                conditions[chr] = {condition:haplotype}
            else:
                conditions[chr][condition] = haplotype

condition_list = list(conditions['3500'].keys())
print(condition_list)

for chr in conditions.keys():
    fig, axs = plt.subplots(3, 1)
    for i in range(0, len(condition_list)):
        hap = conditions[chr][condition_list[i]]
        axs[i].plot(range(len(hap)), hap)
        axs[i].set_title(condition_list[i])

    plt.suptitle('PSB AAGAGGCA-TAGATCGC B-F haplotype sequences ' + chr)
    plt.show()