###################################################################
# Check if the clusters are local or do they overlap on the chromosome?
# Plot a graphs showing the full length of the chromosome and the range
# that each cluster resides in.
###################################################################

import matplotlib.pyplot as plt
import numpy as np
import re

# Initiate variables used for plotting.
count = 0
ranges = []
legend = []

# Read in the cluster data 
with open('../data/conditional_probability_data/clusters_fixed.txt', 'r') as f:
	for line in f:
		count += 1
		cluster_sites = []
		# Find the actual position numbers for each cluster.
		positions = re.findall(r'\d+', line)
		positions = [int(i) for i in positions]
		name = 'Cluster_' + str(count) + '_' + str(len(positions))

		if len(positions) < 10: 	# Ignore small clusters.
			continue

		ranges.append((min(positions), max(positions)))
		legend.append(name)

# Plot each cluster at a separate y axis position.
count = 1
for i in ranges:
	x = range(i[0], i[1], 1)
	y = [count]*len(x)
	plt.plot(x, y, linewidth=5)
	count += 0.5

# Plot the results.
plt.grid()
plt.title('Overlap between big (10+ nodes) clusters (perfect pairs)')
plt.xlabel('chromosome length (bp)')
plt.legend(legend)
plt.show()