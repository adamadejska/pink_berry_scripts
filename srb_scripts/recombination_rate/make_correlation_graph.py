#############################################################################
# Create a graph similar to https://www.nature.com/articles/s41592-018-0293-7
# Where we choose a length l and calculate the probability that we find a mutation
# at distance i+l given that there was a mutation at i
#############################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)
#data = pd.read_csv('/home/ada/Desktop/PinkBerry_scripts_paper/data/srb/SNP_data/srb_snp_data_2024_chr_coverage_6_no_PB93.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

data_index = data.index.values.tolist()[1:]  # names of bacterial samples
chr_list = data.iloc[0, :].tolist()
chr_list = [float(i) for i in chr_list]
data_positions = data.columns.tolist()
data_positions_ints = np.array([int(i.split('.')[0]) for i in data_positions])

all_sigmas = []
# Make a matrix of vectors. Each vector is a binary vector where 1 = different alleles at pos i, and 0 = the same alleles at pos i
for s1i in range(0, len(data_index)):
    for s2i in range(s1i+1, len(data_index)):
        print(s1i, s2i)
        s1, s2 = data_index[s1i], data_index[s2i]
        v1, v2 = np.array(data.loc[s1, :]), np.array(data.loc[s2, :])

        sigma = [s1+'_'+s2]
        for i in range(0, len(v1)):
            
            if v1[i] not in [0, 1] or v2[i] not in [0, 1]:   # not enough info for calculations
                #print(str(v1[i]) + ' ' + str(v2[i]) + ' nan' )
                sigma.append(np.nan)
            elif v1[i] == v2[i]:                     # both alleles are the same
                #print(str(v1[i]) + ' ' + str(v2[i]) + ' same' )
                sigma.append(0)
            else:
                #print(str(v1[i]) + ' ' + str(v2[i]) + ' different' )
                sigma.append(1)                      # the alleles are differen (a mutation and a WT)

        all_sigmas.append(sigma)

tmp_columns = ['name'] + data_positions
sigma_df = pd.DataFrame(all_sigmas, columns=tmp_columns)
sample_names = list(sigma_df.loc[:, 'name'])

# Calculate ds by averaging the correlation matrix
pos_averages = []
for i in range(1, len(tmp_columns)):
    print(i)
    locus_vector = np.array(sigma_df.iloc[:, i])
    numerator = np.nansum(locus_vector)
    denominator = sum(~np.isnan(locus_vector))
    if denominator > 0:
    	pos_averages.append(float(numerator)/float(denominator))

pos_averages = np.array(pos_averages)
ds = sum(pos_averages)/float(len(pos_averages))
ds_var = np.var(pos_averages)

mcorr_out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/mcorr_SRB_data_2024_pipeline_ql_ds_cov_3.csv'
mcorr_out = open(mcorr_out_file, 'w')
mcorr_out.write('# l: the distance between two genomic positions\n')
mcorr_out.write('# m: the mean value of correlatio profile\n')
mcorr_out.write('# v: the variance of correlation profile\n')
mcorr_out.write('# n: the total number of alignments used for calculation\n')
mcorr_out.write('# t: the type of result: Ks is for d_sample, and P2 is for correlation profile\n')
mcorr_out.write('# b: the bootstrap number (all means used all alignments).\n')
mcorr_out.write('l,m,v,n,t,b\n')
mcorr_out.write('0,' + "{:.5f}".format(ds) + ',' + "{:.5f}".format(ds_var) + ',' + str(len(sample_names)) +  ',Ks,all\n')


probabilities_per_distance = {}
for i in range(0, len(data_positions)):
    print(i)
    for j in range(i+1, len(data_positions)):
        pos1, pos2 = data_positions_ints[i], data_positions_ints[j]
        chr1, chr2 = chr_list[i], chr_list[j]

        if abs(int(pos1) - int(pos2)) < 600 and chr1 == chr2:
            dist = abs(int(pos1) - int(pos2))
            pos1_v, pos2_v = list(sigma_df.iloc[:, i]), list(sigma_df.iloc[:, j])

            # Calculate probability of pos2 having 1 given pos1 has 1  
            # P(B|A) = (P and B) / P(A)
            sum_11_positions = 0.0
            sum_1x_positions = 0.0
            total = 0.0

            for k in range(0, len(pos1_v)):
                if pos1_v[k]  in [0, 1] and pos2_v[k] in [0, 1]:
                    total += 1
                    if pos1_v[k] == 1 and pos2_v[k] == 1:
                        sum_11_positions += 1

                    if pos1_v[k] == 1:
                        sum_1x_positions += 1

            if total != 0:
                p_11 = sum_11_positions/total

                if dist not in probabilities_per_distance.keys():
                    probabilities_per_distance[dist] = [p_11]
                else:
                    probabilities_per_distance[dist].append(p_11)


out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/mcorr_SRB_data_2024_ql_ds_cov_3.csv'
out = open(out_file, 'w')

data_out_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/mcorr_SRB_raw_data_2024_ql_ds_cov_3.csv'
raw_out = open(data_out_file, 'w')
raw_out.write('ds,' + str(ds) + '\n')

# Average over all probabilities per distance (sort it as well)
distances_sorted = list(probabilities_per_distance.keys())
distances_sorted.sort()

distances, probabilities = [], []
for k in distances_sorted:
    distances.append(k)
    v = probabilities_per_distance[k]
    raw_out.write(str(k) + ',' + str(v).strip('[').strip(']') + '\n')
    q_l = sum(v)/len(v)
    probabilities.append(q_l / ds)
    out.write(str(k) + ',' + "{:.5f}".format(q_l / ds) + '\n')
    np_v = np.array(v)
    var = np.var(np_v)
    mcorr_out.write(str(k) + ',' + "{:.5f}".format(q_l / ds) + ',' + "{:.5f}".format(var) + ',' + str(len(v)) +',P2,all\n')


out.close()
mcorr_out.close()
raw_out.close()

plt.plot(distances, probabilities, 'o', alpha=0.8, mfc='none')
plt.title('P(l) vs distance SRB (coverage >= 3, no PB93)')
plt.ylabel('average P(l)')
plt.xlabel('Distance l (bp)')
plt.show()
