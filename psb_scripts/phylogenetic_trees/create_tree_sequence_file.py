####################################################################################
# Create files that will be used in the MEGA software for trees creation.
# Take the bacteria SNP data and make a sequence file for each position we are interested in.
# Produces two files: strain names and similarity matrix.
#
# For the PSB data, there is only one possible versions that can be made:
# no 8 clade loci.
####################################################################################

from pysam import VariantFile
import numpy as np
import sys
import pandas as pd



# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)
clade_8_strains = ['PB87','PB40','TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','AAGAGGCA-TATCCTCT','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

# Only consider positions that are not singletons and that are not fixed in our current test clade.
positions = []
counter_used = 0
for i in data_positions:
    
    # Find positions that are singletons -> we don't want those.
    n = sum(np.array(data[i])==1)
    if n > 1:
        positions.append(str(i))
        counter_used += 1

print('used: ' + str(counter_used))
print(len(data_positions))

# Find the indicies for the positions in the dataframe
pos_index = [data_positions.index(str(i)) for i in positions]

strains_kept = []
# Check how many strains have very low number of SNPs
for strain in data_index:
    strain_positions = []
    ns_num = 0.0

    # Change Nans to 0.5 to treat them probabilistically in the distance calculations
    for i in pos_index:
        snp = data.loc[strain].values[i]
        value = 0.5 if np.isnan(snp) else snp
        if value == 0.5:
            ns_num += 1

        strain_positions.append(value)
    
    if float(sum(np.array(strain_positions)==1))/len(strain_positions) >= 0.01:
        strains_kept.append(strain)
    
###
#  Open the bcf file to check which base was at each position of the strain
###

# Read in the original BCF file.
bcf_in = VariantFile('/home/ada/Desktop/Shraiman_lab/data/psb-scaff03.minimap2-freebayes.bcf')  # auto-detect input format
bcf_iter = bcf_in.fetch()
all_samples = list(bcf_in.header.samples)

# Samples that have most loci covered. 
good_samples = ['PB63', 'PB28', 'GGACTCCT-AGAGTAGA', 'AAGAGGCA-AGAGTAGA', 'PB73', 'PB76', 'AAGAGGCA-ACTGCATA', 'GGACTCCT-GTAAGGAG', 'GCTACGCT-GTAAGGAG', 'TCCTGAGC-CTCTCTAT', 'PB55', 'PB_5', 'PB31', 'AAGAGGCA-TATCCTCT', 'PB24', 'PB87', 'AAGAGGCA-AAGGAGTA', 'PB84', 'PB40', 'AGGCAGAA-ACTGCATA', 'AGGCAGAA-CTAAGCCT', 'PB67', 'PB_8', 'PB80', 'PB39', 'PB47', 'PB59', 'PB37', 'PB61', 'TCCTGAGC-TAGATCGC', 'PB77', 'PB69', 'TAAGGCGA-AAGGAGTA', 'PB13', 'PB82', 'PB48', 'PB52', 'GGACTCCT-TAGATCGC', 'PB11', 'GGACTCCT-CTAAGCCT', 'PB44', 'PB36', 'AAGAGGCA-GTAAGGAG', 'PB16', 'PB32', 'PB20', 'GCTACGCT-ACTGCATA', 'PB85', 'PB65', 'PB58', 'PB64', 'CTCTCTAC-GTAAGGAG', 'PB60', 'TAAGGCGA-GTAAGGAG', 'AGGCAGAA-GTAAGGAG', 'PB21', 'TAAGGCGA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'PB78', 'PB_1', 'PB91', 'AAGAGGCA-TAGATCGC', 'CTCTCTAC-AAGGAGTA', 'GTAGAGGA-ACTGCATA', 'CAGAGAGG-AAGGAGTA', 'CTCTCTAC-AGAGTAGA', 'GTAGAGGA-GTAAGGAG', 'PB43', 'PB34', 'AAGAGGCA-CTAAGCCT', 'PB92', 'CTCTCTAC-TATCCTCT', 'AGGCAGAA-AAGGAGTA', 'PB90', 'PB12', 'CAGAGAGG-TATCCTCT', 'PB18', 'PB26', 'GTAGAGGA-AAGGAGTA', 'PB41', 'AGGCAGAA-AGAGTAGA', 'PB75', 'PB45', 'CAGAGAGG-TAGATCGC', 'TAAGGCGA-AGAGTAGA','CAGAGAGG-GTAAGGAG','PB17','CAGAGAGG-ACTGCATA','CTCTCTAC-CTCTCTAT','PB51','CGTACTAG-AGAGTAGA','PB66','PB25','PB83','CGTACTAG-CTAAGCCT','CTCTCTAC-TAGATCGC','TCCTGAGC-AAGGAGTA','CAGAGAGG-AGAGTAGA','PB_3','AGGCAGAA-CTCTCTAT','CAGAGAGG-CTAAGCCT','TAAGGCGA-TAGATCGC','PB81','CGTACTAG-ACTGCATA','CGTACTAG-AAGGAGTA','PB89','TAAGGCGA-ACTGCATA','PB74','CTCTCTAC-ACTGCATA','PB57','PB49','GTAGAGGA-CTAAGCCT','CGTACTAG-CTCTCTAT','GTAGAGGA-CTCTCTAT','PB_9','PB_4','GTAGAGGA-AGAGTAGA','GTAGAGGA-TATCCTCT','PB29','PB10','PB_2','AGGCAGAA-TAGATCGC','CGTACTAG-TATCCTCT','PB53','GGACTCCT-ACTGCATA','PB19','CTCTCTAC-CTAAGCCT','TCCTGAGC-AGAGTAGA','CAGAGAGG-CTCTCTAT','PB35','PB27','GTAGAGGA-TAGATCGC','CGTACTAG-TAGATCGC','TCCTGAGC-CTAAGCCT','TAAGGCGA-CTAAGCCT','TCCTGAGC-TATCCTCT','PB42','TCCTGAGC-GTAAGGAG','TAAGGCGA-CTCTCTAT','TCCTGAGC-ACTGCATA','AGGCAGAA-TATCCTCT','CGTACTAG-GTAAGGAG']

sequences = {}
counter = 0
for pos in bcf_iter:
    #print(pos.pos)
    if str(pos.pos) in positions:
        counter += 1
        print(counter)
        # gather coverage as alleles, mark non-calls, complex alleles and coverage <=3 as 'N'
        tmp_coverage = [-1 if pos.samples[s]['DP'] is None else pos.samples[s]['DP'] for s in good_samples]
        # Get a major allele for each bacteria
        alleles = ['N' if pos.samples[s].alleles[0] is None or len(pos.samples[s].alleles[0])>1 or tmp_coverage[si]<3
                    else pos.samples[s].alleles[0] for si,s in enumerate(good_samples)]
        
        # Add appropriate base to each strains we have kept of the positions we have kept. 
        for s in strains_kept:
            if s not in sequences.keys():
                sequences[s] = alleles[good_samples.index(s)]
            else:
                sequences[s] = sequences[s] + alleles[good_samples.index(s)]



# Save the sequences as fasta file that will be inported to MEGA to create trees.
out = open('/home/ada/Desktop/Shraiman_lab/data/tree_data/PSB_sequences.fasta','w')
for k, v in sequences.items():
    out.write('>' + k + '\n')
    out.write(v + '\n')

out.close()
