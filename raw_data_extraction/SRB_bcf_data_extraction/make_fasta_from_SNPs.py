####################################################
# Make FASTA files for selected strains from of the tree
# based on the whole SNP matrix that will be used 
# in ClonalOrigin analysis. 
####################################################

import numpy as np
import pandas as pd

# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()
data = data.transpose()

# Read in the reference vs alternative allele data
ref_alt_data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/wt_vs_alt_srb.csv')
#print(ref_alt_data)

# Choose 5 strains to represent the different parts of the tree
# TODO: (?) Make a consensus seq for each part of the tree? -> Would probably result in just the consensus.
strains = ['PB50', 'PB_3', 'PB40', 'PB93', 'PB53']       #BasalAB, Clonal, F, Mixing, TopAB
snps = []
dna = ['','','','','']

# Get the SNP data for each chosen strain.
for s in strains:
    snps.append(data[s].tolist())

# For each position in the matrix add the original base or alternative base to make the 
# FASTA sequence.
wt = ref_alt_data['wt'].tolist()
alt = ref_alt_data['alt'].tolist()

for i in range(0, len(snps[0])):
    for s in range(0, len(strains)):
        p = snps[s][i]
        if p == 1:
            dna[s] += alt[i]
        elif p == 0:
            dna[s] += wt[i]
        else:
            dna[s] += 'N'

for s in range(0,len(strains)):
    file_name = '/home/ada/Desktop/Shraiman_lab/srb/clonalorigin_work/data/'+ strains[s] + '_snps_only.fasta'
    with open(file_name, 'w') as f:
        f.write('>'+strains[s]+'\n')
        f.write(dna[s]+'\n')

