###############################################################
# Check the correlation between the distance between pairs of loci and
# the linkage disequilibrium measure (r^2)
# PSB
###############################################################


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# PSB, read the data.
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
data_index = data.index.values.tolist()
data_positions = data.columns.tolist()

# Go through the loci and create pairs. We will calculate the r^2 measure for each pair
# To save time, only concider pairs with a distance between them less than 1 kb. 
out_path = '/home/ada/Desktop/Shraiman_lab/data/LD_pairs_datasets/r_2_PSB_data_2023.csv'
out = open(out_path,'w') 

for i in range(0, len(data_positions)):
    for j in range(i, len(data_positions)):
        if int(data_positions[j]) - int(data_positions[i]) <= 1000:
            locus1 = [str(i).split('.')[0] for i in np.array(data.iloc[:,i])]
            locus2 = [str(i).split('.')[0] for i in np.array(data.iloc[:,j])]

            # Get rid of NaNs
            clean_locus1 = []
            clean_locus2 = []
            for k in range(0, len(locus1)):
                if locus1[k]  in ['0', '1'] and locus2[k] in ['0', '1']:
                    clean_locus1.append(locus1[k])
                    clean_locus2.append(locus2[k])


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
                    out.write(str(data_positions[i]) + ',' + str(data_positions[j]) + ',' + "{:.5f}".format(r_2) + '\n')

out.close()

