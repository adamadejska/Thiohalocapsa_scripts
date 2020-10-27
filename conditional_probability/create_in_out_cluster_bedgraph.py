######################################################################
# Create a bedgraph file to store info about which loci from the 
# original data set are in or outside the big clusters.
# Next, visualize it using IGV program.
######################################################################

import matplotlib.pyplot as plt
import numpy as np
import re

count = 0
in_positions = []

# Read in the cluster data 
with open('../data/conditional_probability_data/conditional_prob_1_clusters.txt', 'r') as f:
	for line in f:
		positions = re.findall(r'\d+', line)
		positions = [int(i) for i in positions]

		if len(positions) < 10: 	# Ignore small clusters.
			continue

		# Add positions to the in-cluster list 
		for i in positions:
			in_positions.append(i)

# Load in the data with all SNP positions
data = np.load('../data/pos_SNPs.npy',allow_pickle='TRUE').item()

out_positions = []
# Look for multi-strain loci that are within clusters'ranges but aren't part of the cluster
inmin = min(in_positions)
inmax = max(in_positions)
for p in data.keys():
	if int(p) not in in_positions and int(p) >= inmin and int(p) <= inmax:
		if sum(np.array(data[p])==1) >= 3:       # Get only multi-strain loci.
			out_positions.append(int(p))
			
# Create an out file to store the data in.
out = open('../data/igv_data/overlap_igv_data_error1_fixed.bedgraph', 'w')

# Create a histogram of SNP number across the chromosome for each cluster
bins = range(inmin, inmax+501, 500)
v, bins, _ = plt.hist(in_positions, bins=bins)

# Write in the positive values in the bedgraph file (inside the cluster)
for j in range(0,len(v)):
	out.write('psb-scaff03' + '\t' + str(bins[j]) + '\t' + str(bins[j+1]) + '\t' + str(v[j]) + '\n')

# Write in the negative values in the bedgraph file (outside the clulster)
v, bins, _ = plt.hist(out_positions, bins=bins)
for j in range(0,len(v)):
	out.write('psb-scaff03' + '\t' + str(bins[j]) + '\t' + str(bins[j+1]) + '\t' + str(v[j]*-1) + '\n')
