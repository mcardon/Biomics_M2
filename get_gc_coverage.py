# -*- coding: utf-8 -*-
# Compute de %GC in high and low coverage regions
# arg1 : ref.fasta (reference genome)
# arg2 : high_cov.csv (csv file of high coverage regions)
# arg3 : low_cov.csv (csv file of low coverage regions)
# arg4 : optional int (minimum length to consider, default = 1)

#import matplotlib.pyplot as plt
import sys
import pandas as pd
from Bio import SeqIO
import collections

filename_ref = str(sys.argv[1])
filename_high = str(sys.argv[2])
filename_low = str(sys.argv[3])
if len(sys.argv) > 4:
	min_len = int(sys.argv[4])
else:
	min_len = 1


# load coverage data
data_hc = pd.read_csv(filename_high)
data_lc = pd.read_csv(filename_low)

# merge into one file
data_hc['coverage'] = pd.Series( ['high' for i in range(data_hc.shape[0])] )
data_lc['coverage'] = pd.Series( ['low' for i in range(data_lc.shape[0])] )
data = data_hc.append(data_lc)

# keep only position col
data = data.loc[:,['chr','start','end','coverage']]
data['length'] = data.loc[:,'end'] - data.loc[:,'start'] +1

# remove duplicates
data = data.drop_duplicates(['chr','start','end'])

# read genome reference
ref_fasta = SeqIO.parse(filename_ref, "fasta")

# get high and low coverage regions
hc_regions = 0
tot_hc = 0
lc_regions = 0
tot_lc = 0


for chromosome in ref_fasta:
	name_chr = chromosome.id
	print("Chromosome : %s" % name_chr)
	data_chr = data.loc[data['chr'] == name_chr]

	for i in range(data_chr.shape[0]):
		# get sequence
		sequence = str(chromosome[ data_chr.iloc[i,1] : (data_chr.iloc[i,2]+1) ].seq)
		tot_len = len(sequence)

		if tot_len > min_len:
			# count nucleotides
			c = collections.Counter(sequence)
			gc_count = c['g'] + c['G'] + c['c'] + c['C']

			if data_chr.iloc[i,3] == 'high' :
				hc_regions += gc_count
				tot_hc += tot_len
			else:
				lc_regions += gc_count
				tot_lc += tot_len

print("For high and low coverage regions of more than %d pb " % min_len)
print("GC percent in high coverage region : %.2f" %(hc_regions/float(tot_hc) ) )
print("GC percent in low coverage region  : %.2f" %(lc_regions/float(tot_lc) ) )
print("GC percent in h/c coverage region  : %.2f" %((hc_regions+lc_regions)/float(tot_hc+tot_lc) ) )










