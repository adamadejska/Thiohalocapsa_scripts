########################################################
# This script looks through a specified gff file and given 
# a list of specified snp positions, it will print 
# which snps are in which genes' ORFs.  
########################################################
import random


# Define a list of SNP positions to look for.
snps = [359318,360162,361305,364122,371925,386455,425329,426665,428251,430617,457513,458149,471595,474996,479317,481530,483528,487388,492072,495979,501174,515443,531511,531819,560742,568507,94501,95825,96827,97057,97112,98440,99413,99432,99474,99515,100326,101252,104128,104350,127277,127312,129965,135115,138285,142442,144513,930509,935677,938478,948055,948401,949689,949914,980584,980792,983338,983910,985642,986242,988422,1478331,1478469,1478494,1478506,1478645,1479407,1480091,1480491,1480525,1480966,1595676,1605010,1620197,1622130,1625512,1627220,1628307,1630075,1630092,2200169,2200346,2200379,2200516,2200841,2201032,2201125,2201215,2201383,2201728,2202241,2202292,2202361,2203476,2212407,2212824,2212834,2214750,2215780,2223748,2223806,2224415,2224763,2300748,2303433,2303858,2304683,2308902,2309283,2313099,2314702,2318754,2320047,2323179,2325537,2342793,2344723,2351085,2389722,2393053,2397080,2397130,2400663,2400970,2401238,2401245,2406112,2406715,2411676,2414666,2417406,2649830,2655413,2661399,2661477,2662655,2664696,2667061,2667723,2682223,2723988,2774899,2780907,2781058,2782156,2793562,2818812,2819789,2830491,2836594,2836816,2836822,2837082,2854573,3121304,3121799,3133725,3133905,3147935,3148854,3151790,3169857,3170448,3170483,3172776,3177239,3178970,3214567,3214776,3235951,3261218,3284124,3284201,3284373,3284382,3290550,3294262,3295390,3297364,3298753,3298877,3987241,3988716,3995821,3997632,3998142,4024954,4026113,4026177,4040052,4045295,4055205,4060379,4062161,4062776,4066945,4069744,4569619,4606607,4626272,4626287,4629459,4633300,4633504,4633509,5046179,5064526,5064625,5064931,5065310,5066064,5066071,5067403,5074441,5075095,5075950,5328865,5329934,5332865,5460184,5461059,5464474,5465736,5465970,5468382,5468502,5468675,5469664,5469936,5470157,5471442,5474448,5475344,5476428,5476570,5476613,5497492,5499031,5500714,5501019,5505125,5505410,6096258,6104551,6104692,6111287,6114402,6119171,6341848,6341883,6355993,6366122,6366161,6369615,6370559,6373602,6438740,6439007,6453004,6457613,6468880,6470796,6476725,6490933,6495911,6499664,6501526,6515552,6515608,7306901,7309780,7310252,7310699,7314004,7361352,7361366,7363913,7381520,7381658,7381718,7539921,7542402,7619363,7826402]
snps.sort()


# Define beginning and end for creating dictionary entries
# That way we save some memory since we're saving a subset of
# the gff file (where the specified snp positions are) rather than 
# the whole thing.
begin = snps[0]-5000
end = snps[-1]+5000
genes = {}
hypothetical_counter = 1
ser_thr_counter = 1
other_counter = 1

# Create a dictionary for gene names and their start and end positions.
with open('../data/igv_files/psb_scaff03_noG.gff', 'r') as f:
	for line in f:
		if line.startswith('#') :
			continue
		elif len(line.split('\t')) != 9:
			continue
		else:
			line = line.split('\t')
			start, stop, description = int(line[3]), int(line[4]), line[-1]
			if start > begin and start < end and stop > begin and stop < end:
				gene_name = description.split(';')[-1].split('=')[-1][:-1]


				# If genes have the same names in multiple places in the genome,
				# create separate entries with name of the gene + a number to
				# make sure the genes aren't overwritten.
				if gene_name == 'hypothetical protein':
					gene_name = gene_name+'_'+str(hypothetical_counter)
					hypothetical_counter += 1
				elif gene_name == 'Serine/threonine-protein kinase pkn1':
					gene_name = gene_name+'_'+str(ser_thr_counter)
					ser_thr_counter += 1
				elif gene_name in genes.keys():
					gene_name = gene_name + '_' + str(other_counter)
					other_counter += 1

				genes[gene_name] = [start, stop]

			if start > end:
				break


# Finally, look through the specified SNPs and check if they are
# inside any of the genes from the dictionary.
# Print the outcome to the console.
for i in snps:
	found = False
	for name in genes:
		if genes[name][0] <= i and genes[name][1] >= i:
			print(str(i) + ' - ' + name)
			found = True
			continue
	if not found:
		print(str(i) + ' - none')