###################################################################################
# Determine clustering in network of highly correlated pairs
# Create a network of highly correlated pairs where distance is over 10 kb
# Look for clusters of pairs in the network to analyze the backbone of the genome.
# We say that the pair is perfect if it follows a perfect vertical descent in the tree.
# joint == min(#_of_i_snps, #_of_j_snps)
###################################################################################

import json
import numpy as np


# Create a graph object using dictionary
graph = {}

# Read the data, determine the distance between pairs, and assign correlation to the dictionary
threshold = 0
# Read in the raw counts for each pair of loci.
with open('../data/conditional_probability_data/conditional_probability_dataset_raw_counts_4plus_fixed.txt', 'r') as f:
	f.readline()   	# Ignore the header
	for line in f:
		i, j, joint, i_snps, j_snps, n = line.split(',')   			# Split the data into each category
		if int(joint) > 0:
			min_val = float(min(int(i_snps), int(j_snps)))
			error = abs(min_val - int(joint))  						# Allow (threshold #) SNPs error margin
			if int(joint) == min(int(i_snps), int(j_snps)):  		# Finding perfectly correlated pairs
			#if error <= threshold:
			
				# Add the connections to the graph (undirected).
				if i not in graph.keys():
					graph[i] = [j]
				else:
					graph[i].append(j)


				if j not in graph.keys():
					graph[j] = [i]
				else:
					graph[j].append(i)

visited = []

#Iterate through the graph to look for clusters using Breadth first search.
for node in graph.keys():
	if node not in visited:
		group = [node]
		to_check = set(graph[node])
		visited.append(node)
		while len(to_check) != 0:
			neighbor = to_check.pop()
			group.append(neighbor) 
			if neighbor not in visited:
				for n in graph[neighbor]:
					if n not in visited:
						to_check.add(n)
			visited.append(neighbor)
			#group.append(neighbor)

# Save the groups into a file.
		print(group)