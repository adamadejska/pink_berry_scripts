####################################################################################
# Create a set of small newick trees with a random set of nodes and certain changes.
# This will be used as an inpuit for Visual TreeCmp to calculate distances between
# the reference tree and the random trees and trees with certain mutations.
# The discance values will be then plotted as a distribution to understand how different
# tree comparison algorithms work and how to interpret the results.
####################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sn

def roundup(n):
    round_num = (n + 9) // 10 * 10
    return(round_num)

# Analysis of the more complex trees.
data_one = pd.read_csv("/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/result_files/VIsual_TreeCmp_simulation_one_cluster_results.csv", index_col=0)
data_two = pd.read_csv("/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/result_files/VIsual_TreeCmp_simulation_two_cluster_results.csv", index_col=0)
data_three = pd.read_csv("/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/result_files/VIsual_TreeCmp_simulation_three_cluster_results.csv", index_col=0)
data_four = pd.read_csv("/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/phylogenetic_trees/tree_comparisons/result_files/VIsual_TreeCmp_simulation_four_cluster_results.csv", index_col=0)

data = [data_one, data_two, data_three, data_four]
methods = ['Triples', 'RFCluster(0.5)', 'MatchingPair', 'NodalSplitted', 'MatchingCluster', 'MAST', 'CopheneticL2Metric']
bin_size = [20, 1, 10, 10, 10, 1, 10]

# For each method, create a distribution containing all of the distances
for i in range(0, len(methods)):
    for j in data:
        method = methods[i]
        distances = np.array(j.loc[:, method])
        total_counts = len(distances)
        max_distance = int(max(distances))
        
        bin_step = bin_size[i]
        bins = range(int(min(distances))-bin_step, roundup(max_distance)+bin_step, bin_step)

        sn.histplot(data=j[method], stat='percent', bins=bins)

    plt.legend(labels=["1 clade moved","2 clades moved", "3 clades moved", "4 clades moved"])
    plt.title(methods[i] + '\n(bin size = ' + str(bin_step) + ' )')
    plt.xlabel('Tree distances')
    plt.ylabel('Percent of random trees')
    plt.show()

