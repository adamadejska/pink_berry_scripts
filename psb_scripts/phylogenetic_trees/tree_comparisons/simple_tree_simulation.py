####################################################################################
# Create a set of small newick trees with a random set of nodes and certain changes.
# This will be used as an inpuit for Visual TreeCmp to calculate distances between
# the reference tree and the random trees and trees with certain mutations.
# The discance values will be then plotted as a distribution to understand how different
# tree comparison algorithms work and how to interpret the results.
####################################################################################

import numpy as np
import random

# Simple tree (reference)
# (A,(B,(C,(D,(E,(F,(G,(H,(I,(J,K))))))))));

"""
node_set = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
n = 1000
random_trees = []

for i in range(0, n):
    random.shuffle(node_set)

    newick_tree = ['(' + node_set[0] + ',(' + node_set[1] + ',(' + node_set[2] + ',(' + node_set[3] + ',(' +
                    node_set[4] + ',(' + node_set[5] + ',(' + node_set[6] + ',(' + node_set[7] + ',(' + 
                    node_set[8] + ',(' + node_set[9] + ',' + node_set[10] + '))))))))));']
    random_trees.append(newick_tree)

for i in random_trees:
    print(i[0])

"""
# More complex tree (reference)
# ((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,(Z,((AA,AB),AC)))))),(B,(((((C,(D,E)),F),G),(H,I)),(J,(K,((L,M),N))))));

# Code for the more complex tree simulations
node_set = ['AD','AE','O','P','R','S','T','U','W','Y','Z','AA','AB','AC','B','C','D','E','F','G','H','I','J','K','L','M','N']
n = 1000
random_trees = []

for i in range(0, n):
    random.shuffle(node_set)

    newick_tree = ['((' + node_set[0] + ',(' + node_set[1] + ',(((((' + node_set[2] + ',(' + node_set[3] + ',' +
                    node_set[4] + ')),' + node_set[5] + '),' + node_set[6] + '),(' + node_set[7] + ',' + 
                    node_set[8] + ')),(' + node_set[9] + ',(' + node_set[10] + ',((' + node_set[11] + ',' + node_set[12]
                    + '),' + node_set[13] + ')))))),(' + node_set[14] + ',(((((' + node_set[15] + ',(' + node_set[16] + ',' + node_set[17]
                    + ')),' + node_set[18] + '),' + node_set[19] + '),(' + node_set[20] + ',' + node_set[21] + ')),(' + node_set[22]
                    + ',(' + node_set[23] + ',((' + node_set[24] + ',' + node_set[25] + '),' + node_set[26] + '))))));']
    random_trees.append(newick_tree)

for i in random_trees:
    print(i[0])

#((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,Z)))),(B,(((((C,(D,E)),F),G),(H,I)),(J,(K,((L,(M,((AA,AB),AC))),N))))));
#((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,(Z,((D,E),C)))))),(B,(((((AC,(AB,AA)),F),G),(H,I)),(J,(K,((L,M),N))))));
#((AD,(AE,(((((O,(P,R)),S),T),(U,W)),(Y,(Z,((D,E),C)))))),(B,(((((N,(L,M)),F),G),(H,I)),(J,(K,((AA,AB),AC))))));

