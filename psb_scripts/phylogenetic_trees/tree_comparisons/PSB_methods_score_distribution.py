####################################################################################
# Create a set of small newick trees with a random set of nodes and certain changes.
# This will be used as an inpuit for Visual TreeCmp to calculate distances between
# the reference tree and the random trees and trees with certain mutations.
# The discance values will be then plotted as a distribution to understand how different
# tree comparison algorithms work and how to interpret the results.
####################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def roundup(n):
    round_num = (n + 9) // 10 * 10
    return(round_num)

# Analysis of the more complex trees.
data = pd.read_csv("/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/PSB_Visual_TreeCmpReport2.csv", index_col=0)

#methods = ['Triples', 'RFCluster(0.5)', 'MatchingPair', 'NodalSplitted', 'MatchingCluster', 'MAST', 'CopheneticL2Metric']
methods = ['MatchingPair', 'NodalSplitted', 'MatchingCluster', 'CopheneticL2Metric']

# For each method, create a distribution containing all of the distances
for i in methods:
    distances = np.array(data.loc[:, i])
    max_distance = int(max(distances))

    if max_distance > 300:
       bin_step = roundup(int(max(np.array(data.loc[:, methods[0]]))//50))
    else:
        bin_step = 1
    
    bin_step = 10
    bins = range(int(min(distances))-bin_step, roundup(max_distance)+bin_step, bin_step)
    plt.hist(distances, bins=bins, edgecolor='black')
    plt.title(i + '\n(bin size = ' + str(bin_step) + ' )')
    plt.xlabel('Tree distances')
    plt.ylabel('Number of random trees')
    plt.show()