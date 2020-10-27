####################################################################
# Check the sizes of the clusters created by the perfect pairs graph.
# Plot them as a histogram.
####################################################################
import matplotlib.pyplot as plt
import numpy as np

# Read in the cluster data.
size = []
with open('../data/conditional_probability_data/clusters_fixed.txt', 'r') as f:
	for line in f:
		line = line.split()
		size.append(len(line))


# Plot the results.
plt.hist(size, bins=range(0,max(size)+1, 1), edgecolor='black')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('size of the cluster (number of nodes)')
plt.ylabel('number of clusters')
plt.title('Cluster sizes for perfect pairs (distance < 50kb) (4+ strains/position)')
plt.show()