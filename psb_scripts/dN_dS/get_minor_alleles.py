# Get minor alleles from the bcf file.

from pysam import VariantFile
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import sys


# Read in the CSV file of the variance data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
good_positions = data.columns.tolist()
good_positions = [int(i) for i in good_positions]
good_samples = index

np.set_printoptions(threshold=sys.maxsize)

bcf_in = VariantFile('/home/ada/Desktop/Shraiman_lab/data/psb-scaff03.minimap2-freebayes.bcf')  # auto-detect input format
bcf_iter = bcf_in.fetch()

# samples that have most loci covered
allele_positions = []
coverage = []
called = []
alleles = []
meta = []
rec = 0
rec_max = 100000
counter = 0
coverage_threshold = 6
for pos in bcf_iter:
	if int(pos.pos)  in good_positions:
		counter += 1
		print(counter)
		# gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
		tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
		# Get a major allele for each bacteria
		tmp_alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<coverage_threshold
					   else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
		tmp_meta = [[tmp_coverage[si], pos.samples[s]['RO'], pos.samples[s]['AO']] for si, s in enumerate(good_samples)]

		# skip if only one proper call or not a proper biallelic
		unique_alleles = set(tmp_alleles)
		if (len(unique_alleles)==2 and 'N' not in unique_alleles) or (len(unique_alleles)==3 and 'N' in unique_alleles):
			alleles.append(tmp_alleles)
			#print(tmp_alleles)
			allele_positions.append(pos.pos)
			coverage.append(tmp_coverage)
			called.append([pos.samples[s].alleles[0] is not None for s in good_samples]) # Not a blank space
			meta.append(tmp_meta)
			rec+=1
			if rec>rec_max:
				break
		else:
			continue
	else:
		continue

# convert to array
alleles = np.array(alleles).T

# make a masked array of SNPs
snps = np.ma.array(alleles, fill_value='N')
snps.mask = (snps=='-')|(snps=='N')  # We want Ns and gaps to stay labeled

# get allele counts and major alleles
alpha = np.array(['A', 'C', 'G', 'T'])
af = np.array([np.sum(snps==a, axis=0) for a in 'ACGT'])

maj = alpha[af.argmax(axis=0)]
print(len(maj))
minor = []

psb_minor = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/dN_dS/psb_minor_alleles_coverage_6.txt'
out_file = open(psb_minor, 'w')
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
                out_file.write(str(allele_positions[i]) + ',' + str(k) + '\n')
                minor.append(k)
    
print('---')
print(len(minor))
print(len(allele_positions))
print(len(good_positions))