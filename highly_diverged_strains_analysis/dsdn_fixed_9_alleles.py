########################################################################
# Synonymous vs nonsynonymous substitution analysis.
# To distinguish between the adaptive and nonadaptive mutations.
# Check if the mutation creates a synonymous or nonsynonymous substitution
# and how does the syn / nonsyn fraction change depending on mutation 
# frequency thresholds.
# Here, the loci positions are specified (dsdn analysis is perfomed on a 
# subset of loci, not the whole file).
########################################################################

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import sys


# Load in important maps and other datasets.

# Dictionary with info about what codon codes for what amino acid.
codon_aa = {'TTT': 'F', 'TTC': 'F',
			'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
			'ATT':'I', 'ATC':'I', 'ATA':'I',
			'ATG':'M',
			'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
			'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
			'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
			'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
			'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
			'TAT':'Y', 'TAC':'Y',
			'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP',
			'CAT':'H', 'CAC':'H',
			'CAA':'Q', 'CAG':'Q',
			'AAT':'N', 'AAC':'N',
			'AAA':'K', 'AAG':'K',
			'GAT':'D', 'GAC':'D',
			'GAA':'E', 'GAG':'E',
			'TGT':'C', 'TGC':'C',
			'TGG':'W',
			'CGT':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
			'AGT':'S', 'AGC':'S',
			'GGT':'G', 'GGA':'G', 'GGC':'G', 'GGG':'G'}

bases = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}  # base pairs used for converting + to - strands.


np.set_printoptions(threshold=sys.maxsize)

# Read in the CSV file of the variance data.
data = pd.read_csv('../data/snp_data_2020.csv', index_col=0)

# Index = names of bacterial samples
index = data.index.values.tolist()
positions = np.array(data.columns.tolist())
nine_fixed = [359318,360162,361305,364122,371925,386455,425329,426665,428251,430617,457513,458149,471595,474996,479317,481530,483528,487388,492072,495979,501174,515443,531511,531819,560742,568507,94501,95825,96827,97057,97112,98440,99413,99432,99474,99515,100326,101252,104128,104350,127277,127312,129965,135115,138285,142442,144513,930509,935677,938478,948055,948401,949689,949914,980584,980792,983338,983910,985642,986242,988422,1478331,1478469,1478494,1478506,1478645,1479407,1480091,1480491,1480525,1480966,1595676,1605010,1620197,1622130,1625512,1627220,1628307,1630075,1630092,2200169,2200346,2200379,2200516,2200841,2201032,2201125,2201215,2201383,2201728,2202241,2202292,2202361,2203476,2212407,2212824,2212834,2214750,2215780,2223748,2223806,2224415,2224763,2300748,2303433,2303858,2304683,2308902,2309283,2313099,2314702,2318754,2320047,2323179,2325537,2342793,2344723,2351085,2389722,2393053,2397080,2397130,2400663,2400970,2401238,2401245,2406112,2406715,2411676,2414666,2417406,2649830,2655413,2661399,2661477,2662655,2664696,2667061,2667723,2682223,2723988,2774899,2780907,2781058,2782156,2793562,2818812,2819789,2830491,2836594,2836816,2836822,2837082,2854573,3121304,3121799,3133725,3133905,3147935,3148854,3151790,3169857,3170448,3170483,3172776,3177239,3178970,3214567,3214776,3235951,3261218,3284124,3284201,3284373,3284382,3290550,3294262,3295390,3297364,3298753,3298877,3987241,3988716,3995821,3997632,3998142,4024954,4026113,4026177,4040052,4045295,4055205,4060379,4062161,4062776,4066945,4069744,4569619,4606607,4626272,4626287,4629459,4633300,4633504,4633509,5046179,5064526,5064625,5064931,5065310,5066064,5066071,5067403,5074441,5075095,5075950,5328865,5329934,5332865,5460184,5461059,5464474,5465736,5465970,5468382,5468502,5468675,5469664,5469936,5470157,5471442,5474448,5475344,5476428,5476570,5476613,5497492,5499031,5500714,5501019,5505125,5505410,6096258,6104551,6104692,6111287,6114402,6119171,6341848,6341883,6355993,6366122,6366161,6369615,6370559,6373602,6438740,6439007,6453004,6457613,6468880,6470796,6476725,6490933,6495911,6499664,6501526,6515552,6515608,7306901,7309780,7310252,7310699,7314004,7361352,7361366,7363913,7381520,7381658,7381718,7539921,7542402,7619363,7826402]
coding_reg = [False]*len(nine_fixed)

pos_mask = []
# Make a map for 9 clade alleles
for i in positions:
	if int(i) in nine_fixed:
		pos_mask.append(True)
	else:
		pos_mask.append(False)

positions = positions[pos_mask]

# Read in the collection sites of the samples.
sites = {}
with open('../data/sample_sites.txt', 'r') as f:
	for line in f:
		sample, site = line.split()
		if sample in index:
			sites[sample] = site

# Get frequency of each SNP position.
all_b = np.zeros(len(positions))
for sample in sites:
	vector = np.array((data.loc[sample].values))
	vector = vector[pos_mask]
	vector[np.isnan(vector)] = 0
	all_b = np.add(all_b, vector)
total_samples = len(sites.keys())



orf_regions = set()
reg_to_descr = {}
# Using the original GFF file, figure out which positions from the SNP matrix hit ORFs.
# Read in the annotation file.
with open('../data/igv_files/psb_scaff03_noG.gff', 'r') as f:
	for line in f:
		# Ignore the header.
		if line.startswith('#') :
			continue
		elif len(line.split('\t')) != 9:
			continue
		else:
			# We found a line with annotated gene. Check if any position from the matrix is inside the gene region.
			line = line.split('\t')
			ann_start, ann_end, strand, descr = int(line[3]), int(line[4]), line[6], line[8].split(';')[-1]
			for i in range(0, len(positions)):
				if int(positions[i]) >= ann_start and int(positions[i]) <= ann_end:
					# Position in the coding region -> we found an ORF.
					coding_reg[i] = True
					print(line)
					if (ann_end+1-ann_start)%3 == 0:
						# Make sure that we are looking at an actually transcribed protein and not an mRNA annotation.
						# (the length will be dividable by 3 (length of a codon))
						orf_regions.add((ann_start, ann_end, strand))
						reg_to_descr[(ann_start, ann_end, strand)] = descr
						break


print('Fraction of coding positions where there are SNPs: %f' %(float(sum(coding_reg))/ len(coding_reg)))
sorted(orf_regions)

# Create a dictionary that will store the data about the ORF start/end positions and the actual FASTA sequence.
orf_to_seq = {}
with open('../data/igv_files/psb_scaff03_fasta_fixed.fasta', 'r') as f:
	for line in f:
		if not line.startswith('>'):
			for region in orf_regions:
				reg_s, reg_e = region[0], region[1]
				orf_to_seq[region] = line[reg_s-1:reg_e]



#  synonymous mutation is a change in the DNA sequence that codes for amino acids 
#  in a protein sequence, but does not change the encoded amino acid.
mutations = {}
for i in positions:
	mutations[int(i)] = ''

# Check what was the mutation (minor allele) for each SNP position.
with open('../data/dsdn/minor_alleles.txt', 'r') as f:
	for line in f:
		tmp = np.array(line.split())
		tmp = tmp[pos_mask]
		for i in range(0, len(tmp)):
			mutations[int(positions[i])] = tmp[i][1:-1]


nonsynonymous = 0
synonymous = 0

# Go through all the sequences to check for mutations
for reg, seq in orf_to_seq.items():
	# Check which SNP positions are contained in the particular sequence
	snp_positions = []
	for i in range(0, len(positions)):
		# if position of SNP is in the region
		if int(positions[i]) >= reg[0] and int(positions[i]) <= reg[1]:
			snp_positions.append(int(positions[i]))

	# Check if we need to reverse the strand. If so, reverse it.
	if reg[2] == '-':
		new_seq = ''
		for i in seq[::-1]:
			new_seq += bases[i]
		seq = new_seq

	# Translate DNA to protein. Look for amino acid with a SNP
	snp_positions = [i-reg[0] for i in snp_positions]
	for i in range(0, len(seq), 3):     # For each codon
		for snp in snp_positions:		# For each snp found in the region
			if snp in range(i, i+3):	# If we're on the good codon

				# Get the original codon and amino acid
				protein = codon_aa[seq[i:i+3]]

				# Get the mutant amino acid
				pos = snp % 3
				mutant = list(seq[i:i+3])
				mutant[pos] = mutations[snp+reg[0]]
				mutant = ''.join(mutant)
				mutant_aa = codon_aa[mutant]

				# Categorize the change as synonymous / nonsynonymous
				if protein != mutant_aa:
					nonsynonymous += 1
					#print('nonsynonymous %d' %(snp+reg[0]))
				elif protein == mutant_aa:
					synonymous += 1
					#print('synonymous %d' %(snp+reg[0]))
				else:
					print('something went wrong')

# Print the information about the raw numbers and also calculated fractions.
print('nonsynonymous: %d' %nonsynonymous)
print('synonymous: %d' %synonymous)
nonsynonymous_f = float(nonsynonymous)/len(positions)
synonymous_f = float(synonymous)/len(positions)
noncoding_f = float(len(positions) - nonsynonymous - synonymous)/len(positions)
if noncoding_f < 0:
	noncoding_f = 0
print('nonsynonymous: %f' %nonsynonymous_f)
print('synonymous: %f' %synonymous_f)
print('nondocing: %f' %noncoding_f)
print('total sites: %d' %len(positions))
