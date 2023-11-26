###############################################################
# Ada Madejska, 2023
# What kind of genes can we find in the SNP-dense hotspots? 
###############################################################

from pysam import VariantFile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in the PSB snp data
gff_path = '/home/ada/Desktop/Shraiman_lab/srb_data/otuB_mp_noG.gff'

counter = 0

# How many genes are the in the gff file?
# Example line:
# unitig_0_otuB|quiver    Prodigal:2.60   CDS     29657   30829   .       +       0       ID=PBSRB_00027;eC_number=2.3.1.16;gene=fadA;inference=ab initio prediction:Prodigal:2.60,similar to AA sequence:UniProtKB:O32177;locus_tag=PBSRB_00027;product=3-ketoacyl-CoA thiolase

with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('unitig'):
            if 'hypothetical protein' in line:
                counter += 1

#print(counter)

# What's the length of the ORFs?
lengths = []
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('unitig'):
            start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
            lengths.append(end-start)

#plt.hist(lengths, bins=range(0,4000, 10))
#plt.show()
"""
# How many contigs are there and how big are they? What's the total length? 
bcf_in = VariantFile('/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf')  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)
counter = 0
current_chr = 'unitig_0_otuB|quiver'
current_pos = 0
all_lengths = []
for pos in bcf_iter:
    tmp_chr = str(pos.chrom)
    if tmp_chr != current_chr:
        print(tmp_chr + ' ' + current_pos)
        all_lengths.append(int(current_pos))
        current_chr = tmp_chr
    else: 
        current_pos = str(pos.pos)

print(tmp_chr + ' ' + current_pos)
print(all_lengths)
all_lengths = np.array(all_lengths)
print(sum(all_lengths))

# How many ORFs per 1 kb? Would expect uniform distribution : 4000 / 46000 -> <1 ORF per 1 kb
# Check the biggest contigs for simplicity sake.
densities = {}
windows = []
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('unitig_0_otuB'):
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
                for i in range(0,l_div+1):
                    int_div = int(start+(i*1000)) // 1000
                    if int_div not in densities.keys():
                        densities[int_div] = 1
                        windows.append(int_div)
                    else:
                        densities[int_div] += 1

plt.hist(densities.values(), bins=range(0,max(densities.values())+1, 1))
plt.show()

"""
# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
contigs_list = [str(i) for i in set(data.iloc[0,:].tolist())]

# [0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
# 20.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 36.0, 37.0, 38.0, 40.0,
# 41.0, 42.0, 43.0, 44.0, 45.0, 49.0, 60.0, 65.0, 66.0, 67.0, 69.0, 70.0, 72.0, 76.0, 77.0, 78.0, 79.0,
# 80.0, 82.0, 83.0, 85.0, 86.0, 87.0, 88.0, 89.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0]

# Count the number of SNPs across all strains in each 1 kb window.
# For now just look at the few biggest contigs: 0, 67, 69, 70, 87
#contigs_list = ['0.0', '67.0', '69.0', '70.0', '87.0']
all_densities = []
all_windows = []
window = 1000
#contigs_list = ['0.0']
mask = [True if x == 'PB93' else False for x in data_index]
for i in contigs_list:
    densities = {}
    windows = []
    for p in range(0,len(data_positions)):
        if str(data.iloc[0,p]) == i:
            #print(data_positions[p])
            # Check if the position has a mutation in PB93 (supermutant) and ignore any singletons from it.
            locus_PB93 = np.array(data[data_positions[p]].tolist())[mask]
            locus_all = np.array(data[data_positions[p]].tolist())
            if sum(locus_PB93==1) == 1 and sum(locus_all==1) == 1:   # check if there's a mutation in PB93 and it's a singleton 
                continue
            else:
                int_div = int(data_positions[p].split('.')[0]) // window
                if int_div not in densities.keys():
                    densities[int_div] = 1
                    windows.append(int_div)
                else:
                    densities[int_div] += 1
    # Look for dense windows ( more than 30 SNPs )
    regions_to_check = []
    for k, v in densities.items():
        if v >= 20:
            regions_to_check.append(int(k))

    # open the gff file and look for genes at those regions
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('unitig_' + str(i.split('.')[0]) + '_'):
                contig = line.split('\t')[0]
                start, end = int(line.split('\t')[3]), int(line.split('\t')[4])
                if int(start // window) in regions_to_check:
                    # retrieve the uniprotkb id
                    info = line.split('\t')[8]
                    info = info.split(';')
                    try:
                        uniprot_index = [i for i, elem in enumerate(info) if 'UniProtKB:' in elem]
                        uniprot_id = info[uniprot_index[0]].split(':')[-1]  # UniProt ID (probably the most important)
                    except IndexError:
                        uniprot_id = 'NA'
                    print( contig + ' ' + str(start) + ' ' + uniprot_id)

    all_densities.append(densities)
    windows = np.array(windows) / 1000
    all_windows.append(windows)



