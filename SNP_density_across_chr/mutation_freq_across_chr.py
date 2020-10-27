###################################################################
# Visualize the biallelic mutations frequency across the genome
# Create a histogram of normalized (SNPs/strain #s) mutations for 
# each geographical group.
# Do that for all the group in all resolutions ( add titles )
# Then do that for all the groups together in all resolutions
###################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import sys

np.set_printoptions(threshold=sys.maxsize)

def create_graph(vector, positions, step, start, end, group, n):
	# Create a dictionary ( position on chr : number of variants )
	
	chrom_d = {}
	for i in range(0, len(vector)):
		if vector[i] != 0:
			real_pos = int(positions[i])
			chrom_d[real_pos] = int(vector[i])

	# NEED TO SPECIFY RESOLUTION
	keys = np.arange(start, end, step)
	sorted_positions = chrom_d.keys()
	sorted_positions.sort()
	clustered = {key: 0 for key in keys}

	# Create a dataframe to store info about the chromosome positions.
	cdf = pd.DataFrame(0.0, index=keys, columns=['n'])
	
	for i in range(0, len(keys)):
		for j in range(0, len(sorted_positions)):
			# Add mutation counts to the right range of positions.
			if sorted_positions[j] >= keys[i] and sorted_positions[j] < keys[i]+step:
					cdf.at[keys[i], 'n'] += chrom_d[sorted_positions[j]]
		cdf.at[keys[i], 'n'] = cdf.at[keys[i], 'n'] / float(n)
		if cdf.at[keys[i], 'n'] >= 6:
			print(keys[i])
	

	# Plot the dataframe as a bar graph.
	cdf.plot(kind='bar', legend=False, width=0.8, figsize=(30,6))
	ticks = cdf.index.values.tolist()
	for i in range(0, len(ticks)):
		if i % 5 != 0:
			ticks[i] = ''
	plt.xticks(range(len(clustered)), ticks, rotation='vertical')
	title = 'Group ' + group + ' biallelic SNPs across chromosome from ' + str(start/1000000) + ' Mb to ' + str(end/1000000) + ' Mb'
	plt.title(title)
	plt.ylabel('Number of SNPs / number of samples')
	plt.xlabel('Position along the chromosome with ' + str(step/1000) + ' Kb resolution')
	plt.show()

# Read in the CSV file of the variance data.
data = pd.read_csv('snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples, positions = chromosomal positions
index = data.index.values.tolist()
positions = data.columns.tolist()

# Read in the collection sites of the samples.
sites = {}
with open('sample_sites.txt', 'r') as f:
	for line in f:
		sample, site = line.split()
		if sample in index:
			sites[sample] = site

# Segregate the samples into their geographical locations.
a, b, c, d, e, f = {}, {}, {}, {}, {}, {}
# Count the number of samples we have in each group.
na, nb, nc, nd, ne, nf = 0, 0, 0, 0, 0, 0
for sample in sites:
	vector = (data.loc[sample].values)
	vector[np.isnan(vector)] = 0
	if sites[sample] == 'A':
		na += 1
		a[sample] = vector
	if sites[sample] == 'B':
		nb += 1
		b[sample] = vector
	if sites[sample] == 'C':
		nc += 1
		c[sample] = vector
	if sites[sample] == 'E':
		ne += 1
		e[sample] = vector
	if sites[sample] == 'F':
		nf += 1
		f[sample] = vector

# Add vectors from one group together
# This will be used to calculate the average SNP amount in each position.
e_vector = np.zeros(len(e[e.keys()[0]]))
for i in e.keys():
	e_vector = np.add(e[i], e_vector)

a_vector = np.zeros(len(a[a.keys()[0]]))
for i in a.keys():
	a_vector = np.add(a[i], a_vector)

b_vector = np.zeros(len(b[b.keys()[0]]))
for i in b.keys():
	b_vector = np.add(b[i], b_vector)

c_vector = np.zeros(len(c[c.keys()[0]]))
for i in c.keys():
	c_vector = np.add(c[i], c_vector)

f_vector = np.zeros(len(f[f.keys()[0]]))
for i in f.keys():
	f_vector = np.add(f[i], f_vector)

graphs = [(a_vector, 'A', na), (b_vector, 'B', nb), (c_vector, 'C', nc), (e_vector, 'E', ne), (f_vector, 'F', nf)]

# Need to average the result over the number of bacteria in a group for good comparison.
# Create graphs in specified regions of the chromosome for each location.
for i in graphs:
	vector = i[0]
	group = i[1]
	n = i[2]
	create_graph(vector, positions, 10000, 0, 2000000, group, n)
	create_graph(vector, positions, 10000, 2000000, 4000000, group, n)
	create_graph(vector, positions, 10000, 4000000, 6000000, group, n)
	create_graph(vector, positions, 10000, 6000000, 8000000, group, n)


"""
############################################
# TEST - looks good
test_v = [1, 6, 5, 0, 0, 3, 2, 8, 9, 1, 0, 1, 3]
t_posi = [2, 15, 19, 21, 23, 31, 33, 41, 50, 58, 67, 68, 79]
create_graph(test_v, t_posi, 10, 0, 90, 'TEST', 2)
"""


