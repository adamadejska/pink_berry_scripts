###############################################################
# Check the correlation between the distance between pairs of loci and
# the linkage disequilibrium measure (r^2)
# SRB
###############################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# SRB, read the data.
# Read in the SRB snp data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021_chr.csv', index_col=0)
data = data.drop(index=('PB93'))    # PB93 is a hypermutator so get rid of it for this analysis

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()
contigs_list = [67, 69, 70, 87]

# Go through the loci and create pairs. We will calculate the r^2 measure for each pair
# To save time, only concider pairs with a distance between them less than 1 kb. 
out_path = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_SRB_data_noPB93_2023_2.csv'
out = open(out_path,'w') 
current_contig = 'x'
for c in contigs_list:
    for i in range(0, len(data_positions)):
        for j in range(i+1, len(data_positions)):
            # Make sure the loci are not super far apart and that they are on the same contig.
            if int(data.iloc[0,i]) == c and int(data.iloc[0,j]) == c:
                if int(data_positions[j].split('.')[0]) - int(data_positions[i].split('.')[0]) <= 1000:
                    
                    if current_contig != int(data.iloc[0,i]):
                        print(int(data.iloc[0,i]))
                        current_contig = int(data.iloc[0,i])


                    locus1 = [str(i).split('.')[0] for i in np.array(data.iloc[1:,i])]
                    locus2 = [str(i).split('.')[0] for i in np.array(data.iloc[1:,j])]

                    # Get rid of NaNs
                    clean_locus1 = []
                    clean_locus2 = []
                    for k in range(0, len(locus1)):
                        if locus1[k]  in ['0', '1'] and locus2[k] in ['0', '1']:
                            clean_locus1.append(locus1[k])
                            clean_locus2.append(locus2[k])

                    if len(clean_locus1) < 20:   # Don't include pairs that don't hold much data. 
                        break

                    # Pair up the alleles
                    alleles = []
                    for k in range(0, len(clean_locus1)):
                        alleles.append(clean_locus1[k]+clean_locus2[k])

                    #print(alleles)
                    # Count the frequency of all allele combinations (00, 01, 10, 11)
                    freq_00 = alleles.count('00')
                    freq_01 = alleles.count('01')
                    freq_10 = alleles.count('10')
                    freq_11 = alleles.count('11')
                    total = float(freq_00 + freq_01 + freq_10 + freq_11)
                    if total != 0:
                        freq_00 = freq_00/total
                        freq_01 = freq_01/total
                        freq_10 = freq_10/total
                        freq_11 = freq_11/total

                        freq_1x = clean_locus1.count('1')
                        freq_0x = clean_locus1.count('0')
                        freq_x1 = clean_locus2.count('1')
                        freq_x0 = clean_locus2.count('0')

                        locus1_total = float(freq_1x + freq_0x)
                        locus2_total = float(freq_x1 + freq_x0)

                        freq_1x = freq_1x / locus1_total
                        freq_0x = freq_0x / locus1_total
                        freq_x1 = freq_x1 / locus2_total
                        freq_x0 = freq_x0 / locus2_total

                        if freq_1x*freq_0x*freq_x1*freq_x0 != 0:
                            r_2 = pow((freq_00*freq_11) - (freq_01*freq_10), 2) / (freq_1x*freq_0x*freq_x1*freq_x0)

                            if r_2 > 1.1:
                                print(r_2)
                                print(alleles)
                                print([freq_00, freq_01, freq_10, freq_11])

                            # Save the data to a file
                            out.write(str(data_positions[i]) + ',' + str(data_positions[j]) + ',' + str(c) + ',' + "{:.5f}".format(r_2) + '\n')

out.close()

