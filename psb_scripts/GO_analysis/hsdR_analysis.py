####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
#
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# 
# Here we focus on XerD since it resolves chromosome dimers formed as 
# a result of DNA replication. There are multiple copies of this genes
# across PSB chromosome.
#
# PSB
####################################################################

from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# How many mutations are there for all XerD? Make a distribution for mutations.
"""
all_mutations = []
gene_code = 'Q7A801'
with open(gff_file_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff03'):
                if gene_code in line:
                    info = line.split('\t')[8]
                    info = info.split(';')

                    start, end = int(line.split('\t')[3]), int(line.split('\t')[4])


                    data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

                    # Index = names of bacterial samples
                    data_index = data.index.values.tolist()
                    data_positions = data.columns.tolist()

                    # Each mutation is an array where each entry represents a specific berry (index).
                    # Make a dictionary for bacteria and count the number of mutations each have per gene.
                    # Inititate the dictionary and set every berry to 0. 
                    berry_specific_mutations = {}
                    for berry in data_index:
                            berry_specific_mutations[berry] = 0

                    mutation_counter = 0
                    for i in data_positions:
                        if int(i) >= start and int(i) <= end:
                            mutation_counter += 1
                    all_mutations.append(mutation_counter)

plt.hist(all_mutations, edgecolor='black')
plt.xlabel('number of mutations per gene copy')
plt.ylabel('number of gene copies')
plt.title('Distribution of number of mutations per hsdR gene copy')
plt.show()
"""

gene_code = 'P60240'
with open(gff_file_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff03'):
                if gene_code in line:
                    info = line.split('\t')[8]
                    info = info.split(';')

                    start, end = int(line.split('\t')[3]), int(line.split('\t')[4])


                    data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

                    # Index = names of bacterial samples
                    data_index = data.index.values.tolist()
                    data_positions = data.columns.tolist()

                    # Each mutation is an array where each entry represents a specific berry (index).
                    # Make a dictionary for bacteria and count the number of mutations each have per gene.
                    # Inititate the dictionary and set every berry to 0. 
                    berry_specific_mutations = {}
                    for berry in data_index:
                            berry_specific_mutations[berry] = 0

                    mutation_counter = 0
                    for i in data_positions:
                        if int(i) >= start and int(i) <= end:
                            mutation_counter += 1
                            tmp = np.array(data[i])
                            for j in range(0, len(data_index)):
                                if not np.isnan(tmp[j]):
                                    berry_specific_mutations[data_index[j]] += tmp[j]

                    # Take only the berries that have mutations. We don't want to plot a bunch of zeros.
                    valid_berries, num_of_mutations = [], []
                    for k,v in berry_specific_mutations.items():
                        if v > 0:
                            valid_berries.append(k)
                            num_of_mutations.append(v)

                    num_of_mutations, valid_berries = zip(*sorted(zip(num_of_mutations, valid_berries)))

                    plt.barh(valid_berries, num_of_mutations)
                    plt.title('Type-1 restriction enzyme R protein hsdR ' + str(start) + ' - ' + str(end) + '\n number of mutations: ' + str(mutation_counter))
                    plt.xlabel('Number of mutations per berry')
                    plt.ylabel('Berry names')
                    plt.grid(axis = 'x')
                    plt.show()
"""

# What is the number of mutations for each berry IN TOTAL. Meaning, if we look at all copies of the XerD
# gene together, how many mutations would each berry have?
gene_code = 'Q7A801'

data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Inititate the dictionary and set every berry to 0. 
berry_specific_mutations = {}
for berry in data_index:
        berry_specific_mutations[berry] = 0


with open(gff_file_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff03'):
                if gene_code in line:
                    info = line.split('\t')[8]
                    info = info.split(';')

                    start, end = int(line.split('\t')[3]), int(line.split('\t')[4])

                    # Each mutation is an array where each entry represents a specific berry (index).
                    # Make a dictionary for bacteria and count the number of mutations each have per gene.

                    mutation_counter = 0
                    for i in data_positions:
                        if int(i) >= start and int(i) <= end:
                            mutation_counter += 1
                            tmp = np.array(data[i])
                            for j in range(0, len(data_index)):
                                if not np.isnan(tmp[j]):
                                    berry_specific_mutations[data_index[j]] += tmp[j]


# Take only the berries that have mutations. We don't want to plot a bunch of zeros.
valid_berries, num_of_mutations = [], []
for k,v in berry_specific_mutations.items():
    if v > 1:    # Don't look at sinlge mutations in a berry to make the graph more readable
        valid_berries.append(k)
        num_of_mutations.append(v)

num_of_mutations, valid_berries = zip(*sorted(zip(num_of_mutations, valid_berries)))

plt.barh(valid_berries, num_of_mutations)
plt.title('Type-1 restriction enzyme R protein hsdR \n all copies together' )
plt.xlabel('Number of mutations per berry')
plt.ylabel('Berry names')
plt.grid(axis = 'x')
plt.show()
"""