# -*- coding: utf-8 -*-
# arg1 : reference fasta file
# arg2 : bed file
# arg3 : int : window size min (optionnal : default 101)
# arg4 : int : window size max (optionnal : default 103)
# arg5 : int : step (optionnal : default 10)

import matplotlib.pyplot as plt
import sys
import numpy as np
from sequana import *


filename_ref = str(sys.argv[1])
filename_bed = str(sys.argv[2])

# filename_ref = "/home/mcardon/Mel/Datas/Ref_genome_NC/NC_002929_Bpertussis.fasta"
# filename_bed = "/home/mcardon/Mel/Datas/Ref_genome_NC/cov_pacbio_blasr_mapq100.bed"


if len(sys.argv) > 3:

	# if even number, add 1 to get odd number
	window_min = int(sys.argv[3])
	if not window_min % 2:
		window_min += 1

	window_max = int(sys.argv[4])

	# if odd number, add 1 to get even number
	step = int(sys.argv[5])
	if step % 2:
		step += 1

else:
	window_min = 101
	window_max = 103
	step = 10


cov = GenomeCov(filename_bed)
# compute gc content with default window

for window in range(window_min, window_max, step):
	print("Compute GC content with window %d" %window)
	cov.compute_gc_content(filename_ref,circular = True, window_size=window)

	for chromosome in cov.chr_list:
		coverage_data = chromosome.df.loc[:,'cov']
		mu = np.mean(coverage_data)
		sigma = np.std(coverage_data)

		# plot all data
		corr = chromosome.plot_gc_vs_coverage()
		plt.title("%s - Window = %d - Cor = %.3f" % (chromosome.chrom_name,window, corr))
		plt.savefig('Cor_GC_coverage_window_%d_%s.png' % (window, chromosome.chrom_name ))
		#plt.show()
		plt.close()

		# plot 95% data
		corr = chromosome.plot_gc_vs_coverage()
		plt.title("%s - Window = %d - Cor = %.3f (6 std interval)" % (chromosome.chrom_name,window, corr))
		plt.xlim(mu - 6*sigma, mu + 6*sigma)
		plt.savefig('Cor_GC_coverage_window_%d_%s_6std.png' % (window, chromosome.chrom_name ))
		#plt.show()
		plt.close()

