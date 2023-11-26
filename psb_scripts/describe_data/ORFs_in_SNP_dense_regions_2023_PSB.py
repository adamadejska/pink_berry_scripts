###############################################################
# Ada Madejska, 2023
# What kind of genes can we find in the SNP-dense hotspots? 
###############################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in the PSB snp data
gff_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

counter = 0

# How many genes are the in the gff file?
# Example line:
# psb-scaff03     Prodigal:2.6    CDS     77728   78417   .       -       0       ID=P1s_00083;eC_number=3.1.3.18;Name=gph_1;gene=gph_1;inference=ab initio prediction:Prodigal:2.6,protein motif:HAMAP:MF_00495;locus_tag=P1s_00083;product=Phosphoglycolate phosphatase
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff'):
            if 'hypothetical protein' in line:
                counter += 1

print(counter)

# What's the length of the ORFs?
lengths = []
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff'):
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
            lengths.append(end-start)

plt.hist(lengths, bins=range(0,5000, 10))
plt.show()

# How many ORFs per 1 kb? Would expect uniform distribution : 7000 / 8000 -> <1 ORF per 1 kb
densities = {}
windows = []
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff'):
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
            length = end-start
            if length < 1000:
                int_div = int(start) // 1000
                if int_div not in densities.keys():
                    densities[int_div] = 1
                    windows.append(int_div)
                else:
                    densities[int_div] += 1

            else:
                l_div = length // 1000
                for i in range(0,l_div):
                    int_div = int(start+(i*1000)) // 1000
                    if int_div not in densities.keys():
                        densities[int_div] = 1
                        windows.append(int_div)
                    else:
                        densities[int_div] += 1

plt.hist(densities.values(), bins=range(0,max(densities.values())+1, 1))
plt.show()

###
# Do the same analysis for the PSB 
### 

data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
gff_path = '/home/ada/Desktop/Shraiman_lab/data/psb_scaff03_noG.gff'

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Count the number of SNPs across all strains in each 1 kb window.
window = 1000
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

# Keep track on which regions are SNP dense.
dense_windows = []
for k,v in densities.items():
    if v > 10:
        dense_windows.append(k)

print(dense_windows)

# Check those dense windows against the gff file to check what kind of genes are there.
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('psb-scaff') and 'hypothetical protein' not in line:
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
            if start // window in dense_windows:
                # retrieve the uniprotkb id
                info = line.split('\t')[8]
                info = info.split(';')
                try:
                    uniprot_index = [i for i, elem in enumerate(info) if 'UniProtKB:' in elem]
                    uniprot_id = info[uniprot_index[0]].split(':')[-1]  # UniProt ID (probably the most important)
                except IndexError:
                    uniprot_id = 'NA'
                print( str(start) + ' ' + uniprot_id)

