###########################################################################
# Empirical cumulative distribution of distances between loci fixed in 9 clade.
# This script checks how the fixed loci in the 9 clade are distributed in other genomes.
# There are two options: They are either evenly distributed along 
# the chromosome or they are clustered.
# By checking the distances between them we can check which is the case.
###########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random


def find_min_distances(snps):
# This function takes a list of genomic snp positions.
# For each snp it looks at its immediate neighbours and stores the min distances from the two.
# It returns the list of minimal distances between neighbouring snps.

	min_distances = []

	for i in range(1,len(snps)-1):
		# Calculate the distances.
		dist1 = abs(int(snps[i])-int(snps[i+1]))
		dist2 = abs(int(snps[i])-int(snps[i-1]))

		# Choose minimal distance and store it.
		d = min([dist1,dist2])
		min_distances.append(d)

	return(min_distances)


def make_controls(snp_num, data_positions):
# This function makes random set of snps from all snp data that has the same number as genome we're looking at
# Then, for each random snp set, it calculates the minimal distances.
# Those random controls will be used to see if the minimal distances we find from the genomes are significantly different.

	control_num = 10
	control_distances = []  # list of lists of controls

	for i in range(0,control_num):
		# Choose x randomly sampled snps from the dataset with not repeats. 
		random_snps = random.sample(list(data_positions),snp_num)
		random_snps.sort()

		# Calculate the minimal distance between neighbouring snps
		min_distances = find_min_distances(random_snps)

		# Save the min distances for plotting.
		control_distances.append(min_distances)

	return(control_distances)


def make_plot_variables(min_dist):
# This function creates the x and y variables for the plot function.
# The x axis is the minimal distance we see in a particular sample. 
# The y axis is the relative amount of data that we've seen below a certain x value.

	n = 1.0/len(min_dist)  # length of each vertical step. (Y axis goes from 0 to 1)
	y, x = [], []   # vectors to store y and x values
	min_dist_set = list(set(min_dist))
	min_dist_set.sort()

	for i in range(0,len(min_dist_set)):
		x.append(min_dist_set[i])

		if i == 0:
			y_previous = 0
		else:
			y_previous = y[i-1]

		y_val = y_previous + len(np.where(np.array(min_dist)==min_dist_set[i])[0])*n
		y.append(y_val)

	return(x,y)

# Main
# Positions that are fixed in the 9 clade. Will be used for filtering loci in other strains.
nine_clade_fixed_loci=[47637,47684,48680,48708,52690,52884,65806,75861,75983,78801,78807,94501,95825,96827,97057,97112,98440,99413,99432,99474,99515,100326,101252,104128,104350,127277,127312,129965,135115,138285,142442,144513,145422,169507,169638,172490,180534,185193,187656,190768,198334,227235,249288,268589,318695,331162,331694,331782,333319,337509,345634,359318,360162,361305,364122,371925,386455,425329,426665,428251,430617,457513,458149,471595,474996,479317,481530,483528,487388,492072,495979,501174,515443,531511,531819,560742,568507,590154,650973,652354,731030,733130,733685,734017,736352,739600,746242,746263,749252,752793,758683,766559,767771,767780,768584,776574,780396,780933,780997,781048,781075,783992,784045,794585,794824,825861,859110,859480,859885,859985,862327,864048,864132,869201,872328,872542,877557,897700,930509,935677,938478,948055,948401,949689,949914,980584,980792,983338,983910,985642,986242,988422,996796,1001291,1008468,1012061,1012656,1024963,1027262,1034512,1034732,1037830,1039364,1040899,1051375,1054452,1054639,1055488,1055494,1059461,1059606,1062597,1062615,1065312,1065877,1070291,1070466,1072319,1093425,1133962,1134396,1134613,1137018,1143991,1156064,1157216,1162894,1168602,1173498,1177851,1178643,1178660,1187073,1189323,1191341,1192649,1193825,1195476,1207899,1252999,1253009,1270003,1270061,1271442,1271457,1271838,1271847,1271865,1271882,1271892,1272318,1272648,1274255,1274280,1274342,1274769,1274874,1274926,1277395,1288411,1289284,1291926,1292108,1301367,1309506,1315554,1321219,1321333,1323781,1326168,1329282,1340522,1365024,1372910,1373185,1373191,1388341,1391645,1392067,1392072,1392095,1392525,1392539,1396209,1399929,1449655,1453715,1459755,1478331,1478469,1478494,1478506,1478645,1479407,1480091,1480491,1480525,1480966,1518431,1523784,1535511,1537600,1547596,1548361,1558164,1565291,1574557,1583555,1584398,1595676,1605010,1620197,1622130,1625512,1627220,1628307,1630075,1630092,1650750,1659260,1669837,1679148,1708533,1708539,1714946,1716488,1716733,1719750,1728416,1736289,1736907,1742717,1742723,1765722,1765920,1766686,1769046,1778699,1779274,1779329,1779677,1779873,1780609,1780766,1781183,1781192,1783283,1783544,1794848,1794872,1794899,1795319,1797680,1797987,1799537,1800132,1801802,1802014,1802128,1804912,1805031,1805589,1805670,1805976,1806283,1818781,1823994,1825564,1828513,1848157,1848548,1862121,1862803,1898950,1902803,1909152,1909594,1923188,1930033,1931905,1953583,1976448,1978929,2003831,2004997,2007519,2018834,2019308,2022298,2046611,2051090,2057362,2058812,2059235,2079563,2084457,2086510,2086566,2087192,2089389,2092696,2129220,2163459,2164873,2171922,2173725,2173791,2181656,2186097,2190354,2190747,2191032,2191572,2192280,2193402,2193690,2193996,2194252,2200169,2200346,2200379,2200516,2200841,2201032,2201125,2201215,2201383,2201728,2202241,2202292,2202361,2203476,2212407,2212824,2212834,2214750,2215780,2223748,2223806,2224415,2224763,2229922,2231430,2242316,2243522,2244088,2245740,2247819,2249273,2300748,2303433,2303858,2304683,2308902,2309283,2313099,2314702,2318754,2320047,2323179,2325537,2342793,2344723,2351085,2389722,2393053,2397080,2397130,2400663,2400970,2401238,2401245,2406112,2406715,2411676,2414666,2417406,2432152,2432234,2448146,2458513,2471185,2471550,2473081,2477219,2480353,2480671,2484010,2493838,2608783,2611868,2613601,2614293,2614572,2615286,2615873,2616688,2616713,2617635,2649830,2655413,2661399,2661477,2662655,2664696,2667061,2667723,2682223,2723988,2741038,2751528,2754944,2755051,2760938,2764011,2764257,2773469,2774899,2780907,2781058,2782156,2793562,2797682,2804809,2804887,2805021,2808258,2808294,2808513,2811302,2811347,2813113,2818812,2819789,2830491,2836594,2836816,2836822,2837082,2854573,2864770,2872306,2921855,2922267,2935290,2964016,2965598,2966368,2972658,2979727,2997071,3009923,3027563,3031777,3032950,3032956,3034187,3036785,3037145,3039214,3039567,3039588,3039810,3039828,3039876,3040026,3040035,3040113,3040129,3040224,3040595,3040692,3040802,3040942,3041048,3041104,3041150,3041167,3041231,3041315,3041426,3041651,3041687,3041693,3041843,3044463,3047387,3052357,3053852,3077094,3081377,3082605,3084231,3086166,3099432,3107069,3107211,3108664,3112593,3115433,3121304,3121799,3133725,3133905,3147935,3148854,3151790,3169857,3170448,3170483,3172776,3177239,3178970,3191328,3198401,3199100,3201968,3214567,3214776,3235951,3261218,3284124,3284201,3284373,3284382,3290550,3294262,3295390,3297364,3298753,3298877,3432463,3433224,3481332,3481482,3482006,3484352,3484368,3484870,3485299,3485380,3485680,3489916,3498981,3502290,3520596,3603732,3603753,3627324,3632557,3634385,3634777,3635276,3644764,3645048,3651611,3683468,3683700,3686176,3702771,3702957,3703998,3708835,3708973,3709090,3709261,3709270,3709366,3710556,3749235,3756460,3764313,3767109,3767116,3768107,3777009,3781320,3806292,3808323,3812034,3818950,3825418,3831577,3834364,3851207,3865057,3884444,3884460,3884532,3887365,3891535,3895944,3896018,3898602,3901243,3911038,3921077,3977788,3984510,3987241,3988716,3995821,3997632,3998142,4024954,4026113,4026177,4040052,4045295,4055205,4060379,4062161,4062776,4066945,4069744,4125760,4128575,4129022,4136514,4141210,4141223,4142201,4143055,4187860,4188148,4305495,4435180,4436036,4436054,4456453,4457275,4460609,4460719,4464759,4471966,4475825,4477552,4477625,4493248,4496366,4501974,4507191,4508445,4515463,4516356,4521715,4525616,4527217,4528845,4560883,4569619,4606607,4626272,4626287,4629459,4633300,4633504,4633509,4635517,4674185,4674293,4686973,4716853,4718284,4718353,4724741,4740987,4740996,4751987,4753728,4754987,4757699,4757774,4760506,4774172,4781809,4781815,4784849,4786019,4786184,4787077,4789655,4794301,4794314,4794501,4794525,4796144,4797298,4801638,4803776,4804374,4808485,4809442,4809511,4809576,4813615,4816141,4816473,4817678,4818608,4819175,4819805,4824653,4827805,4829720,4836063,4836600,4838731,4860676,4861382,4863156,4863220,4867429,4871368,4880363,4885157,4889002,4899583,4905252,4906898,4906947,4907282,4908611,4909919,4910432,4910572,4911088,4911906,4914371,4918735,4923168,4923196,4936477,4954449,4956901,4991208,4991294,4991623,4992578,4992621,4994424,4994449,4994940,4995068,4995081,4995121,4995149,4995222,4996771,4996984,4998497,4998712,4998870,4998875,5046179,5064526,5064625,5064931,5065310,5066064,5066071,5067403,5074441,5075095,5075950,5093083,5105215,5105945,5113539,5170193,5177243,5204390,5211148,5211429,5215932,5223884,5234541,5238266,5241889,5242783,5243131,5249142,5260569,5299072,5299503,5299513,5300080,5328865,5329934,5332865,5346067,5355402,5370601,5372576,5373565,5377232,5378458,5379874,5382064,5382092,5393817,5424719,5460184,5461059,5464474,5465736,5465970,5468382,5468502,5468675,5469664,5469936,5470157,5471442,5474448,5475344,5476428,5476570,5476613,5497492,5499031,5500714,5501019,5505125,5505410,5510395,5515686,5515749,5517029,5534604,5552504,5558677,5562797,5588654,5588975,5589236,5589398,5589648,5596012,5596450,5596573,5596899,5599753,5604191,5604950,5606107,5612217,5641352,5641463,5643654,5644037,5644091,5645508,5651793,5655124,5655142,5655149,5655273,5655285,5660822,5662276,5674536,5680025,5680780,5680987,5682240,5683894,5683924,5684497,5684710,5684723,5690348,5690403,5690944,5691951,5693030,5693739,5694794,5696406,5795891,5802644,5818452,5944990,5966308,5985873,6012209,6012241,6012569,6012626,6013103,6013145,6013255,6013273,6013283,6013306,6074751,6091729,6096258,6104551,6104692,6111287,6114402,6119171,6153964,6154850,6155880,6157025,6157276,6161073,6165257,6169424,6172315,6172402,6176862,6189152,6189922,6190730,6190801,6190848,6198543,6206397,6236787,6238374,6239091,6242708,6244164,6250589,6255724,6272563,6341848,6341883,6355993,6366122,6366161,6369615,6370559,6373602,6438740,6439007,6453004,6457613,6468880,6470796,6476725,6490933,6495911,6499664,6501526,6515552,6515608,6518022,6541607,6541683,6544711,6547964,6571778,6571820,6573850,6575539,6584643,6586998,6588713,6593428,6630523,6630586,6630610,6630616,6631446,6631747,6631864,6631954,6631975,6637583,6638129,6642211,6642361,6642579,6642992,6643829,6649934,6650299,6720756,6722942,6739112,6748060,6754044,6770367,6771661,6778381,6778974,6779881,6780477,6787980,6793863,6794027,6794273,6798915,6801929,6802594,6802731,6802862,6802929,6803189,6803195,6803228,6803299,6803472,6808607,6836883,6838771,6891646,6906109,6937234,6940778,6945225,6949537,6969841,6971786,7015806,7017403,7018608,7018936,7022527,7022532,7022649,7022706,7025380,7039306,7043721,7062193,7096377,7097227,7112866,7113038,7113253,7114961,7141458,7144796,7149773,7152524,7153260,7161752,7163840,7165477,7167616,7168828,7169363,7171592,7171746,7172008,7172584,7172618,7172652,7172690,7172747,7173014,7173078,7173777,7173815,7173935,7173998,7174037,7174043,7174171,7174184,7174193,7174199,7174211,7174364,7174371,7174376,7174681,7176194,7176200,7176222,7176292,7176459,7176477,7176528,7176946,7177899,7177908,7178005,7179003,7181791,7182512,7182531,7182676,7182774,7182811,7183163,7183317,7183326,7183555,7183797,7184149,7184389,7184600,7184728,7184859,7185056,7186873,7201703,7204509,7205159,7206369,7206672,7206751,7206761,7206876,7206957,7207617,7211159,7211175,7211657,7211834,7212050,7306901,7309780,7310252,7310699,7314004,7315392,7316362,7319458,7320710,7322121,7323248,7325340,7328163,7332339,7333411,7333590,7336038,7337352,7346502,7353319,7361352,7361366,7363913,7381520,7381658,7381718,7404954,7405155,7407837,7414136,7418765,7421252,7421624,7421813,7422707,7469575,7473947,7476954,7477314,7477356,7477425,7477461,7477592,7477598,7477606,7512634,7539921,7542402,7619363,7826402,7941423,7941429,7941456,7947437]

# Load the original data matrix.
og_data = pd.read_csv('../data/snp_data_2020.csv', index_col=0)
og_index = og_data.index.values.tolist()   # list of strain names
data_positions = og_data.columns.tolist()  # list of genome position indices
data_positions = np.array([int(i) for i in data_positions])

# List of strains that belong to the highly diverged 9 clade.
nine_clade = ['TCCTGAGC-CTCTCTAT','TCCTGAGC-TAGATCGC','PB73','PB87','PB40','GGACTCCT-CTAAGCCT','GGACTCCT-GTAAGGAG','GGACTCCT-AGAGTAGA','PB80']

# Go through all the strains.
for strain in og_index:
	# Get all SNP data for particular strain
	strain_snps = np.array(list(og_data.loc[strain]))

	# Filter the SNP data by mutation. If strain has a mutation at pos x, keep that; discard everything else.
	mutations_pos = data_positions[list(np.where(strain_snps == 1)[0])]

	# Filter mutated positions by 9 clade fixed loci. We are only interested in the positions that are also fixed in the 9 clade.
	fixed_mut_pos = [i for i in mutations_pos if i in nine_clade_fixed_loci]

	# Sort mutated positions just in case.
	fixed_mut_pos.sort()

	# Some strains don't have any fixed 9 clade loci. Ignore those strains.
	if len(fixed_mut_pos) <= 20:
		continue

	# Calculate minimal distances between neighbouring SNPs
	min_distances = find_min_distances(fixed_mut_pos)

	# Make controls that will be used to compare the genome positions to see if they are significant.
	controls = make_controls(len(fixed_mut_pos), nine_clade_fixed_loci)

	# Plot the genome distance and control distances as a empirical cumulative distribution.
	strain_x, strain_y = make_plot_variables(min_distances)

	# Calculate average control cumulant.
	all_controls = []
	for c in controls:
		all_controls += c

	control_x,control_y = make_plot_variables(all_controls)
		

	
	# Find the max distance between two distribiutions to test for significance.
	s,c = 0,0
	max_vert_dist, max_x_val = 0,0
	for i in range(0,strain_x[-1]):
		# Check each step in the distribution 
		if i in strain_x:
			# If we see step in the genome distance vallues, then we need to update the y value for comparisons
			s = strain_y[strain_x.index(i)]

		if i in control_x:
			c = control_y[control_x.index(i)]

		if abs(s-c) > max_vert_dist:
			max_vert_dist = abs(s-c)
			max_x_val = i

	print('max vertical distance: ' + str(max_vert_dist))
	print('max vertical distance at genomic distance: ' + str(max_x_val))

	# Make a plot of strain and averaged controls
	plt.plot(control_x,control_y,'b', label='averaged control')
	plt.plot(strain_x,strain_y, 'r', label='strain data')
	plt.axvline(x=max_x_val, color='green', alpha=0.5)

	title = 'Minimal distances between neighbouring SNPs cumulative distribution \n' + strain + ' (' + str(len(fixed_mut_pos)) + ' fixed loci)'
	plt.title(title)
	plt.ylabel('Fraction of distances')
	plt.xlabel('Distance (bp)')
	plt.legend()

	# Save the figures to a file
	fig_name = 'strains_cumulants_dist_btwn_fixed_loci/cumulant_' + strain + '.png'
	plt.savefig(fig_name)
	plt.clf()
	#plt.show()