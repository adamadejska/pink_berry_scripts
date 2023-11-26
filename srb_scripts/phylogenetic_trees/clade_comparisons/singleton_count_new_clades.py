####################################################################################
# After creating a new tree without HGTed loci, we want to check how "old" are the 
# newly created clades. 
# 
# Do the clades that contain strains from the mixing layer have more singletons than the rest?
# SRB analysis.
####################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_boxplots(df):
    """
    Makes a boxplot based on a pandas dataframe.
    df: pandas dataframe with rows as samples and columns as different conditions 
    """
    vals, names, xs = [],[],[]
    for i, col in enumerate(df.columns):
        vals.append(df[col].values)
        names.append(col)
        xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0])) 
        # adds jitter to the data points - can be adjusted

    plt.boxplot(vals, labels=names)
    palette = ['magenta', 'orange', 'green', 'blue', 'red', 'teal']

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)


# Read in the SNP data
data = pd.read_csv('/home/ada/Desktop/Shraiman_lab/srb_data/srb_snp_data_2021.csv', index_col=0)

# Delete PB93 because it's a hyper mutator
data = data.drop('PB93')

clade1 = ['GGACTCCT-TAGATCGC', 'GGACTCCT-ACTGCATA', 'PB11', 'CTCTCTAC-TAGATCGC', 'PB77', 'PB35', 'PB29', 'GGACTCCT-CTAAGCCT', 'CAGAGAGG-CTCTCTAT', 'CTCTCTAC-CTCTCTAT', 'TCCTGAGC-CTAAGCCT', 'PB85', 'CAGAGAGG-CTAAGCCT']
clade2 = ['PB40', 'PB88', 'PB16', 'PB64', 'AAGAGGCA-TAGATCGC', 'AAGAGGCA-GTAAGGAG', 'PB87']
clade3 = ['CAGAGAGG-AAGGAGTA', 'CAGAGAGG-GTAAGGAG', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'AAGAGGCA-CTAAGCCT', 'GCTACGCT-ACTGCATA', 'PB39', 'CAGAGAGG-TAGATCGC', 'CAGAGAGG-AGAGTAGA', 'CAGAGAGG-ACTGCATA']
clade4 = ['CTCTCTAC-CTAAGCCT', 'CTCTCTAC-ACTGCATA', 'PB37', 'PB_5', 'CAGAGAGG-TATCCTCT', 'PB45', 'CTCTCTAC-AGAGTAGA']
clade5 = ['AGGCAGAA-AGAGTAGA', 'PB66', 'PB43', 'CGTACTAG-ACTGCATA', 'TCCTGAGC-TATCCTCT', 'CGTACTAG-TAGATCGC', 'TAAGGCGA-AAGGAGTA', 'CGTACTAG-AGAGTAGA', 'CGTACTAG-CTCTCTAT', 'PB_3', 'PB44', 'PB41', 'TCCTGAGC-GTAAGGAG', 'PB25']
clade6 = ['TAAGGCGA-ACTGCATA', 'TAAGGCGA-CTCTCTAT', 'AGGCAGAA-CTCTCTAT', 'AGGCAGAA-TAGATCGC', 'PB83', 'PB33', 'GTAGAGGA-CTCTCTAT', 'PB53']

whole_tree_no_hgt = ['GGACTCCT-TAGATCGC', 'GGACTCCT-ACTGCATA', 'PB11', 'CTCTCTAC-TAGATCGC', 'PB77', 'PB35', 'PB29', 'GGACTCCT-CTAAGCCT', 'CAGAGAGG-CTCTCTAT', 'CTCTCTAC-CTCTCTAT', 'TCCTGAGC-CTAAGCCT', 'PB85', 'CAGAGAGG-CTAAGCCT', 'PB40', 'PB88', 'PB16', 'PB64', 'AAGAGGCA-TAGATCGC', 'AAGAGGCA-GTAAGGAG', 'PB87', 'CAGAGAGG-AAGGAGTA', 'CAGAGAGG-GTAAGGAG', 'CTCTCTAC-TATCCTCT', 'AAGAGGCA-TATCCTCT', 'AAGAGGCA-CTCTCTAT', 'AAGAGGCA-CTAAGCCT', 'GCTACGCT-ACTGCATA', 'PB39', 'CAGAGAGG-TAGATCGC', 'CAGAGAGG-AGAGTAGA', 'CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT', 'CTCTCTAC-ACTGCATA', 'PB37', 'PB_5', 'CAGAGAGG-TATCCTCT', 'PB45', 'CTCTCTAC-AGAGTAGA', 'PB57', 'PB21','AGGCAGAA-AGAGTAGA', 'PB66', 'PB43', 'CGTACTAG-ACTGCATA', 'TCCTGAGC-TATCCTCT', 'CGTACTAG-TAGATCGC', 'TAAGGCGA-AAGGAGTA', 'CGTACTAG-AGAGTAGA', 'CGTACTAG-CTCTCTAT', 'PB_3', 'PB44', 'PB41', 'TCCTGAGC-GTAAGGAG', 'PB25', 'PB_4', 'TAAGGCGA-ACTGCATA', 'TAAGGCGA-CTCTCTAT', 'AGGCAGAA-CTCTCTAT', 'AGGCAGAA-TAGATCGC', 'PB83', 'PB33', 'GTAGAGGA-CTCTCTAT', 'PB53','TCCTGAGC-AGAGTAGA', 'PB84', 'PB58', 'TCCTGAGC-ACTGCATA', 'PB42', 'CGTACTAG-CTAAGCCT']

data_index = data.index.values.tolist()  # names of bacterial samples
data_positions = data.columns.tolist()

list_of_clades = [clade1, clade2, clade3, clade4, clade5, clade6]
clade_names = ['Clade 1', 'Clade 2', 'Clade 3', 'Clade 4', 'Clade 5', 'Clade 6']
clade_singleton_dict = {}

for j in range(0, len(list_of_clades)):

    clade = list_of_clades[j]
    singletons_in_clade = [np.NaN] * 15    # to save the number of singletons in each strain of the clade

    print('starting clade ' + clade_names[j])
    print(len(clade))

    for k in range(0, len(clade)):
        
        strain = clade[k]
        singletons = 0          # reset the counter

        for i in data_positions:
        
            # Find positions that are singletons -> we want to count those.
            n = sum(np.array(data[i])==1)
            if n == 1 and data.loc[strain, i] == 1:
                singletons += 1
        
        singletons_in_clade[k] = singletons

    clade_singleton_dict[clade_names[j]] = singletons_in_clade
    print(singletons_in_clade)

    print('done with clade ' + clade_names[j])

df = pd.DataFrame(clade_singleton_dict)      
plot_boxplots(df)
plt.yscale('log')
plt.ylabel('Number of singletons')
plt.title('Number of singletons per strain in each SRB clade')
plt.show()
