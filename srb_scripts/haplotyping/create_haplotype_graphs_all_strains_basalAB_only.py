###########################################################################
# Haplotyping strains from the mixing layer.
# This script attempts to haplotype the strains from the mixing layer based 
# on three haplotype profiles: F clade (G), basal_AB clade (Y), and wild type (WT).
# First, it creates the probabilistic profiles for each haplotype and next it
# scans a given strain using a window and calculating a probability that 
# the particular window comes from a particular haplotype.
###########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Code from: https://www.geeksforgeeks.org/bellman-ford-algorithm-dp-23/
# Python3 program for Bellman-Ford's single source
# shortest path algorithm.
# Class to represent a graph
class Graph:
 
    def __init__(self, vertices):
        self.V = vertices # No. of vertices
        self.graph = []
 
    # function to add an edge to graph
    def addEdge(self, u, v, w):
        self.graph.append([u, v, w])
         
    # utility function used to print the solution
    def printArr(self, dist):
        print("Vertex Distance from Source")
        for i in range(self.V):
            print("{0}\t\t{1}".format(i, dist[i]))
     
    # The main function that finds shortest distances from src to
    # all other vertices using Bellman-Ford algorithm. The function
    # also detects negative weight cycle
    def BellmanFord(self, src):
 
        # Step 1: Initialize distances from src to all other vertices
        # as INFINITE
        dist = [float("Inf")] * self.V
        dist[src] = 0
        p = [-1] * self.V      # Solve what nodes create the shortest path
 
 
        # Step 2: Relax all edges |V| - 1 times. A simple shortest
        # path from src to any other vertex can have at-most |V| - 1
        # edges
        for _ in range(self.V - 1):
            # Update dist value and parent index of the adjacent vertices of
            # the picked vertex. Consider only those vertices which are still in
            # queue
            for u, v, w in self.graph:
                if dist[u] != float("Inf") and dist[u] + w < dist[v]:
                        dist[v] = dist[u] + w
                        p[v] = u
                         
        # Print shortest distance to the last node
        #print('shortest distance to last node: %d' %dist[self.V-1])
        
        
        # Print path from start to the last node
        path = []
        v = self.V-1
        while v != -1:
            path.append(v)
            v = p[v]

        path.reverse()
        return(path)
 
 
def calculate_haplotype_probability(snp_matrix):
    """
    This function calculates a probability of a haplotype given a snp matrix.
    It counts the number of each allele and divides it by all the alleles.
    The haplotypes are: green, yellow, WT, missing.
    """
    haplotype_counts = [0, 0, 0]  # G, W, M
    strains = snp_matrix.index.tolist()

    for s in strains:
        seq = np.array(snp_matrix.loc[s,:])
        for i in seq:
            if i == 1:
                haplotype_counts[0] += 1
            elif i == -1:
                haplotype_counts[1] += 1
            elif i == 0:
                haplotype_counts[2] += 1
    
    total = float(len(strains)) * len(snp_matrix.columns.tolist())
    haplotype_fractions = np.array(haplotype_counts) / total
    
    return(haplotype_fractions)


def make_profile_matrix(strains, snp_matrix):
    """
    This funxtion creates a probability matrix that's 4 x m where m is the number of columns 
    (SNP positions) in the data matrix and 4 indicates the alphabet 
    (green allele (1), yellow allele (2), white allele (-1), missing allele (0))
    The output is a matrix filled with probabilities of each allele being found at each position.
    """

    df = snp_matrix.loc[strains,:]
    columns = df.columns.tolist()
    prob_df = pd.DataFrame(index=['G','W','M'], columns=columns)

    # Calculate the probability of seeing each allele at each position in the matrix.
    # Add pseudocounts so that we don't end up with probability of 0.
    for i in columns:
        col = np.array(df.loc[:,i])
        g = (sum(col==1)+1)/(float(len(strains))+3)
        w = (sum(col==-1)+1)/(float(len(strains))+3)
        m = (sum(col==0)+1)/(float(len(strains))+3)

        prob_df[i]['G'] = g
        prob_df[i]['W'] = w
        prob_df[i]['M']= m

    return(prob_df)

def make_type_profile_matrix(strains, snp_matrix, type):
    """
    This function makes a matrix where most of the probability is on white allele.
    We add pseudocounts to get rid of the possibility of log(0) calculations.
    It returns a matrix filled with probabilities for each position in the SNP matrix.
    """
    df = snp_matrix.loc[strains,:]
    columns = df.columns.tolist()
    prob_df = pd.DataFrame(index=['G','W','M'], columns=columns)

    # Calculate the probability of seeing each allele at each position in the matrix.
    # Add pseudocounts so that we don't end up with probability of 0.
    for i in columns:
        g = (1)/(float(len(strains))+3)
        if type == 'WT':
            w = (len(strains))/(float(len(strains))+3)
            m = (1)/(float(len(strains))+3)
        elif type == 'M':
            w = (1)/(float(len(strains))+3)
            m = (len(strains))/(float(len(strains))+3)

        prob_df[i]['G'] = g
        prob_df[i]['W'] = w
        prob_df[i]['M']= m

    return(prob_df)


def calculate_window_probability(window, g_profile, w_profile, m_profile):
    """
    This function calculates the probability of the sequence of the window based on its location
    in the sequence and the haplotype profiles.
    P(sequence | Haplotype) = product(P(allele i | H))  for i alleles in the window
    """
    sym_to_hap = {-1:'W', 0:'M', 1:'G'} 
    hap_to_index = {'G':0,'W':1,'M':2}

    # Calculate the likelihood for the whole window
    g_prob, w_prob, m_prob = 1.0, 1.0, 1.0
    for j in range(0, len(window)):
        g_prob *= g_profile.iloc[hap_to_index[sym_to_hap[window[j]]], j]
        w_prob *= w_profile.iloc[hap_to_index[sym_to_hap[window[j]]], j]
        m_prob *= m_profile.iloc[hap_to_index[sym_to_hap[window[j]]], j]

    return([g_prob, w_prob, m_prob])


def calculate_sequence_probability(haplotype_fractions, g_profile, w_profile, m_profile, window):
    """
    This function calculates the probability of a sequence for each haplotype to discriminate the 
    regions of no discrimination. 
    P(s) = sum(P(s|H)P(H))   for all haplotypes
    haplotype_fractions = [G, Y, W, M]
    """
    window_probabilities = calculate_window_probability(window, g_profile, w_profile, m_profile)
    total_probability = 0.0
    for i in range(0, len(haplotype_fractions)):
        total_probability += window_probabilities[i] * haplotype_fractions[i]

    return(total_probability)


def calculate_penalty(columns, max_penalty):
    """
    This function calculates the penalty of switching between haplotypes based on the 
    real genomic distance between SNPs. The longer the distance, the switching becomes more 
    probable. 
    It returns a list of penalties based on distances.
    """
    penalties = []

    for i in range(0, len(columns)-1):
        dist = abs((float(float(columns[i+1])) - int(float(columns[i]))) / 1000) # make it into kb

        # Calculate penalty based on power law function penatly = 8*e^(-x/L)
        l = 7
        p = max_penalty*np.exp(-1/l*dist)
        if p > max_penalty:   # Cap it
            p = max_penalty

        penalties.append(p)
    
    return(penalties)


def make_all_profiles(basal_AB_strains, snp_matrix):

    # Make profile probability matrices for each haplotype
    green_prob_matrix = make_profile_matrix(basal_AB_strains, snp_matrix)
    white_prob_matrix = make_type_profile_matrix(basal_AB_strains, snp_matrix, 'WT')  # strain input doesn't really matter here.
    missing_prob_matrix = make_type_profile_matrix(basal_AB_strains, snp_matrix, 'M')

    return([green_prob_matrix, white_prob_matrix, missing_prob_matrix])


def run_haplotyping(start, length, window_size, strain, contig, snp_matrix, max_penalty, profiles):
    # Info on haplotypes and dataset
    basal_AB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']

    keep_green = True
    if strain in basal_AB_strains:
        keep_green = False

    # Take only the columns that belong to the contig
    columns = snp_matrix.columns.tolist()[start:length]

    haplotype_fractions = calculate_haplotype_probability(snp_matrix) # G, W, M

    green_prob_matrix, white_prob_matrix, missing_prob_matrix = profiles


    #test_strain = 'CAGAGAGG-AGAGTAGA'
    seq = np.array(snp_matrix.loc[strain,:])

    # Scan through the sequence using a window (size n=10 for now but up to optimization).
    # For each window calculate the probability using the profiles that this window came from a particular
    # profile. Choose the profile with highest probability

    green_prob = []
    white_prob = []
    missing_prob = []
    for i in range(start, length - window_size):
        window = np.array(seq[i:i+window_size])
        g_profile = green_prob_matrix.iloc[:,i:i+window_size]
        w_profile = white_prob_matrix.iloc[:,i:i+window_size]
        m_profile = missing_prob_matrix.iloc[:,i:i+window_size]

        seq_prob = calculate_sequence_probability(haplotype_fractions, g_profile, w_profile, m_profile, window)
        g, w, m = calculate_window_probability(window, g_profile,  w_profile, m_profile)

        if keep_green:
            green_prob.append(np.log((g*haplotype_fractions[0])/seq_prob))
        else:
            green_prob.append(-10)

        white_prob.append(np.log((w*haplotype_fractions[1])/seq_prob))
        missing_prob.append(np.log((m*haplotype_fractions[2])/seq_prob))



    plt.clf()
    # Make a log likelihood graph.
    plt.plot(range(length-window_size-start),green_prob, c='green', label='F haplotype')
    plt.plot(range(length-window_size-start),white_prob, c='gray', label='WT haplotype', alpha=0.6)
    plt.plot(range(length-window_size-start),missing_prob, c='lightgreen', label='missing haplotype')

    plt.xlabel('ordered SNP sequence (bp)')
    plt.ylabel('log likelihood')

    plt.legend()
    plt.title(strain + ' contig ' + str(contig) + ' log likelihood graph for haplotypes')
    plt.tight_layout()
    plt.grid()
    #plt.show()

    ########
    # Use the Bellman-Ford algorithm to create a shortest path using dynamic programming. 
    # The graph is used to choose which haplotype will create a path with the smallest penalty score
    # at the end. This will represent our most probable shifts from one haplotype to another.
    ######## 

    g = Graph(2+(len(green_prob)*3))
    # Add start node (start, end, weight)
    g.addEdge(0, 1, green_prob[0])
    g.addEdge(0, 2, white_prob[0])
    g.addEdge(0, 3, missing_prob[0])

    contig = int(contig.split('.')[0])
    penalties = calculate_penalty(columns, max_penalty)
    green_nodes, white_nodes, missing_nodes = [], [], []

    # add all other edges and nodes
    for i in range(0,len(green_prob)-1):

        green_node_current = (3*i)+1
        white_node_current = (3*i)+2
        missing_node_current = (3*i)+3

        green_nodes.append(green_node_current)
        white_nodes.append(white_node_current)
        missing_nodes.append(missing_node_current)

        green_node_next = (3*(i+1))+1
        white_node_next = (3*(i+1))+2
        missing_node_next = (3*(i+1))+3

        # Add a green node and all of the edges from it. 
        # (start, end, weight)
        g.addEdge(green_node_current, green_node_next, abs(green_prob[i+1]))
        #print('g.addEdge(%d, %d, %f)' %(green_node_current, green_node_next, abs(green_prob[i+1])))
        g.addEdge(green_node_current, white_node_next, penalties[i])
        #print('g.addEdge(%d, %d, %f)' %(green_node_current, white_node_next, penalty))
        g.addEdge(green_node_current, missing_node_next, penalties[i])

        # Add a white node and all of the edges form it.
        g.addEdge(white_node_current, white_node_next, abs(white_prob[i+1]))
        #print('g.addEdge(%d, %d, %f)' %(white_node_current, white_node_next, abs(white_prob[i+1])))
        g.addEdge(white_node_current, green_node_next, penalties[i])
        #print('g.addEdge(%d, %d, %f)' %(white_node_current, green_node_next, penalty))
        g.addEdge(white_node_current, missing_node_next, penalties[i])

        # Add a missing node and all of the edges from it.
        g.addEdge(missing_node_current, missing_node_next, abs(missing_prob[i+1]))
        #print('g.addEdge(%d, %d, %f)' %(yellow_node_current, yellow_node_next, abs(yellow_prob[i+1])))
        g.addEdge(missing_node_current, white_node_next, penalties[i])
        #print('g.addEdge(%d, %d, %f)' %(yellow_node_current, white_node_next, penalty))
        g.addEdge(missing_node_current, green_node_next, penalties[i])
        

    # Add a sink node
    #print(1+(len(green_prob)*3))
    g.addEdge(green_node_next, 1+(len(green_prob)*3), green_prob[len(green_prob)-1])
    g.addEdge(white_node_next, 1+(len(green_prob)*3), white_prob[len(white_prob)-1])
    g.addEdge(missing_node_next, 1+(len(green_prob)*3), missing_prob[len(missing_prob)-1])

    # Calculate the solution starting from the source node.
    path = g.BellmanFord(0)

    # Translate the solution into a graph.
    haplotype = []
    for i in path:
        if i in green_nodes:
            haplotype.append(1)
        elif i in white_nodes:
            haplotype.append(0)
        elif i in missing_nodes:
            haplotype.append(-2)

    plt.clf()
    plt.plot(range(len(haplotype)), haplotype)
    plt.title(strain + ' contig ' + str(contig) + ' B-F haplotype sequence')
    #plt.show()
    return(haplotype)


def get_haplotype_lengths(haplotype, positions):
    """
    Given a list of consecutive SNP matrix positions and the haplotype we figured out each position belongs to,
    figure out the length of each region for the green and yellow haplotype. 
    """

    green_lengths = []

    # Green haplotype
    haplotype += [0]
    green_hap_only = np.array(haplotype)==1
    green_hap_only = green_hap_only.astype(int)

    region = []
    for i in range(0, len(green_hap_only)):
        if green_hap_only[i] == 1:
            region.append(i)
        elif len(region):
            start = str(positions[region[0]]).split('.')[0]
            end = str(positions[region[-1]]).split('.')[0]
            #print(start, end)
            green_lengths.append(abs(int(end)-int(start)))
            region = []

    return(green_lengths)


if __name__ == "__main__":

    # Read in the SNP matrix and figure out which columns correspond to which contig.
    #snp_matrix_path = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_basalAB_ONLY_threshold50_NaNs30wholegenome_fixed.csv'
    snp_matrix_path = '/home/ada/Desktop/Shraiman_lab/srb_data/snp_matrix_fixed_F_ONLY_threshold50_NaNs30wholegenome_fixed.csv'
    basal_AB_strains = ['CGTACTAG-TATCCTCT','CGTACTAG-GTAAGGAG','PB50','AGGCAGAA-TATCCTCT','GTAGAGGA-TAGATCGC','TAAGGCGA-CTAAGCCT','PB82','PB90','PB27','AGGCAGAA-GTAAGGAG','CGTACTAG-CTAAGCCT','TAAGGCGA-TATCCTCT','PB18','CGTACTAG-AAGGAGTA','GTAGAGGA-AGAGTAGA','PB10','PB34','PB76']
    f_strains = ['CAGAGAGG-AAGGAGTA','CAGAGAGG-GTAAGGAG','CTCTCTAC-TATCCTCT','AAGAGGCA-CTCTCTAT','AAGAGGCA-TATCCTCT','AAGAGGCA-CTAAGCCT','GCTACGCT-ACTGCATA','PB39','PB32','PB40','PB88','PB16','PB79','PB87','AAGAGGCA-TAGATCGC','AAGAGGCA-GTAAGGAG','PB64']

    mixing_layer_strains = ['CAGAGAGG-TAGATCGC','CAGAGAGG-AGAGTAGA','CAGAGAGG-ACTGCATA','CTCTCTAC-CTAAGCCT','CTCTCTAC-ACTGCATA','PB37','PB_5','CAGAGAGG-TATCCTCT','CTCTCTAC-AGAGTAGA','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','PB93','PB11','CTCTCTAC-TAGATCGC','GGACTCCT-TAGATCGC','GGACTCCT-ACTGCATA','PB77','PB29','PB35','PB45','TCCTGAGC-CTAAGCCT','PB85','CAGAGAGG-CTAAGCCT','GGACTCCT-CTAAGCCT','PB57','PB21','TAAGGCGA-CTCTCTAT']
    #mixing_layer_strains = ['PB93']

    backwards = False
    max_penalty = 5

    # Read in the data matrix
    snp_matrix = pd.read_csv(snp_matrix_path, index_col=0)
    columns = snp_matrix.columns.tolist()
    strains = snp_matrix.index.tolist()[1:]

    if backwards:
        # Reverse the snp matrix
        snp_matrix = snp_matrix.transpose()
        snp_matrix = snp_matrix[::-1].reset_index(drop=True)
        snp_matrix = snp_matrix.transpose()
        columns = columns[::-1]
        snp_matrix.columns = columns


    #profiles = make_all_profiles(basal_AB_strains, snp_matrix)
    profiles = make_all_profiles(f_strains, snp_matrix)

    contigs = list(snp_matrix.iloc[0])[:-1]  # for F clade cuz last col is from a different contig
    contig_start_stop = {}

    current_contig = contigs[0]
    for i in range(0, len(contigs)):
        if str(contigs[i]) not in contig_start_stop:
            contig_start_stop[str(contigs[i])] = [i]  # start location of the contig
            if contigs[i] != current_contig:
                contig_start_stop[str(contigs[i-1])].append(i-1)  # stop location of the contig
                current_contig = contigs[i]

    contig_start_stop[str(contigs[i-1])].append(i)  # stop location of the contig

    all_green_lengths = []
    per_strain_blocks_green = {}
    
    #for strain in strains:
    for strain in mixing_layer_strains:   # testing
        all_contigs_green = []
        for contig, vals in contig_start_stop.items():
            if vals[1] - vals[0] > 100 :
                print(strain + ' ' + str(contig))
                #print(contig_start_stop[contig][0],contig_start_stop[contig][1])
                #if str(contig) == '70.0':
                haplotype = run_haplotyping(contig_start_stop[contig][0],contig_start_stop[contig][1], 3, strain, contig, snp_matrix, max_penalty, profiles)

                green_lengths = get_haplotype_lengths(haplotype,columns[contig_start_stop[contig][0]: contig_start_stop[contig][1]])

                all_green_lengths += green_lengths
                
                all_contigs_green += green_lengths

        per_strain_blocks_green[strain] = all_contigs_green

    #print(all_green_lengths)
    #print(all_yellow_lengths)
    print(per_strain_blocks_green)
