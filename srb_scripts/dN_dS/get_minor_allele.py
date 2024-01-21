# Get minor alleles from the bcf file.
# SRB version

from pysam import VariantFile
import pandas as pd
import numpy as np 

srb_bcf = '/home/ada/Desktop/Shraiman_lab/srb_data/srb.minimap2-freebayes.bcf'
srb_data = '/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_3_no_PB93.csv'

# Read in the CSV file of the variance data.
data = pd.read_csv(srb_data, index_col=0)

index = data.index.values.tolist()[1:]     # Index = names of bacterial samples
positions = data.columns.tolist()
positions = [int(float(i)) for i in positions]
contigs = data.iloc[0].tolist()

# Read in the bcf file for the mutation information
bcf_in = VariantFile(srb_bcf)  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

alleles = []
allele_contigs = []
allele_positions = []

counter = 0
for pos in bcf_iter:
    chr = pos.chrom.split('_')[1]
    for i in range(0, len(positions)):
        if int(pos.pos) == positions[i] and int(chr) == contigs[i]:
            counter += 1
            print(counter)
            # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
            tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in index]
            # Get a major allele for each bacteria
            tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<3
                        else pos.samples[s].alleles[0] for si,s in enumerate(index)]

            # skip if only one proper call or not a proper biallelic
            unique_alleles = set(tmp_alleles)
            alleles.append(tmp_alleles)
            allele_contigs.append(chr)
            allele_positions.append(pos.pos)
            
        else:
            continue

print('---')
# convert to array
alleles = np.array(alleles).T

# make a masked array of SNPs
snps = np.ma.array(alleles, fill_value='N')
snps.mask = (snps=='-')|(snps=='N')  # We want Ns and gaps to stay labeled

# get allele counts and major alleles
alpha = np.array(['A', 'C', 'G', 'T'])
af = np.array([np.sum(snps==a, axis=0) for a in 'ACGT'])

maj = alpha[af.argmax(axis=0)]
minor = []

srb_minor = '/home/ada/Desktop/PinkBerry_scripts_paper/srb_scripts/dN_dS/srb_minor_alleles_coverage_3.txt'
out_file = open(srb_minor, 'w')
# Figure out the counts of each allele
for i in range(0, len(af[0])):
    d = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    d['A'] = af[0][i]
    d['C'] = af[1][i]
    d['G'] = af[2][i]
    d['T'] = af[3][i]
    #print(d)
    max_int = d[max(d, key=d.get)]
    for k,v in d.items():
        if k != maj[i]:
            if v <= max_int and v > 0:
                out_file.write(str(allele_positions[i]) + ',' + str(allele_contigs[i]) + ',' + str(k) + '\n')
                minor.append(k)
    
print('---')
print(len(minor))
print(len(allele_contigs))
print(len(allele_positions))
"""
srb_minor = '/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/main_data_files/srb_minor_alleles.txt'
srb_major = '/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/main_data_files/srb_major_alleles.txt'

# Save the minor and major alleles to separate files. (each position comma-separated)
out_file = open(srb_minor, 'w')
minor = ','.join(minor)
out_file.write(minor)
out_file.close()

out_file = open(srb_major, 'w')
major = ','.join(maj)
out_file.write(major)
out_file.close()
"""