#############################################################################
# This script makes a bar graph where x axis is a list of strains in order 
# that they appear on the tree and y axis is the number of fixed 9 clade loci
# in each strain on y axis.
#############################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Define the order of strains and other variables.
tree_order = ['TCCTGAGC-TAGATCGC','TCCTGAGC-CTCTCTAT','PB73','GGACTCCT-GTAAGGAG','GGACTCCT-AGAGTAGA','GGACTCCT-CTAAGCCT','PB40','PB87','AAGAGGCA-TATCCTCT','PB80','PB63','AAGAGGCA-TAGATCGC','PB39','PB8','GCTACGCT-GTAAGGAG','PB24','PB55','AAGAGGCA-AAGGAGTA','AAGAGGCA-ACTGCATA','AAGAGGCA-AGAGTAGA','PB31','PB64','GCTACGCT-ACTGCATA','AAGAGGCA-GTAAGGAG','PB16','PB32','AAGAGGCA-CTAAGCCT','AAGAGGCA-CTCTCTAT','PB47','PB48','PB77','PB45','PB61','PB78','GGACTCCT-ACTGCATA','GGACTCCT-TAGATCGC','PB69','TAAGGCGA-TATCCTCT','PB13','PB85','PB37','PB5','AGGCAGAA-CTAAGCCT','AGGCAGAA-ACTGCATA','PB28','PB67','PB76','PB59','PB52','TAAGGCGA-AAGGAGTA','PB84','AGGCAGAA-AGAGTAGA','PB44','PB36','PB58','PB26','PB34','CAGAGAGG-CTAAGCCT','CTCTCTAC-CTAAGCCT','CAGAGAGG-GTAAGGAG','CAGAGAGG-CTCTCTAT','CTCTCTAC-CTCTCTAT','GTAGAGGA-CTAAGCCT','CTCTCTAC-ACTGCATA','CTCTCTAC-TAGATCGC','CTCTCTAC-TATCCTCT','CAGAGAGG-AGAGTAGA','CAGAGAGG-AAGGAGTA','CAGAGAGG-ACTGCATA','CAGAGAGG-TAGATCGC','CTCTCTAC-AAGGAGTA','CTCTCTAC-AGAGTAGA','CAGAGAGG-TATCCTCT','CTCTCTAC-GTAAGGAG','PB11','PB60','PB90']
fixed_loci_number = []
data_dict = {}


with open('Cluster_46_3107_names_9clade.csv', 'r') as f:
	# There's just one line with names and values separated by commas
	# The name consits of two parts (strainname_numberoffixed9cladeloci) which we want to extract.
	for line in f:
		data = line.split(',')
		for i in data:
			name,val = i.split('_')
			data_dict[name] = int(val)


for i in tree_order:
	fixed_loci_number.append(data_dict[i])


# Plot the results as a bar graph
plt.bar(tree_order, fixed_loci_number)
plt.title('Number of fixed 9 clade SNPs in the strains (based on order in the tree)')
plt.ylabel('Amount of fixed 9 clade SNPs present')
plt.xlabel('Names of strains ordered by the tree')
plt.xticks(rotation=90)
plt.show()