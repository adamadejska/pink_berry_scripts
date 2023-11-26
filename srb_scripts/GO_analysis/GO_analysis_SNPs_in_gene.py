####################################################################
# Use the background set of GO terms from all genes in the gff file
# to check if experimental set are overrepresented in certain terms. 
#
# How many SNPs are in certain genes of interest? 
# Given a gene code, find how many SNPs and where they are. 
# PSB
####################################################################

from collections import Counter
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

gene_code = ['Q8ZNY9', 'D0ZRB2', 'Q0P9C7']   # bio processes PSB
gene_code = ['P26394', 'D0ZRB2', 'P18143', 'P51961', 'Q8FDQ2', 'P12995', 'Q8ZNY9']  # molecular functions PSB
gene_code = 'Q8ZNY9'

gff_file_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# Example gff entry
# psb-scaff03	Prodigal:2.6	CDS	1660	2769	.	+	0	ID=P1s_00002;Name=dnaN;dbxref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q9I7C4;locus_tag=P1s_00002;product=Beta sliding clamp

# Specific to PSB data
with open(gff_file_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff03') and gene_code in line:
            info = line.split('\t')[8]
            info = info.split(';')

            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])

print(start, end)
print("length: " + str(end-start))

data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

num_of_mutations = 0
num_of_singletons = 0
for i in data_positions:
    if int(i) >= start and int(i) <= end:
        num_of_mutations += 1
        tmp = np.array(data[i])
        print(np.array(data[i]))
        print(np.nansum(tmp))
        if np.nansum(tmp) == 1:
            num_of_singletons += 1

print(num_of_mutations)
print(num_of_singletons)