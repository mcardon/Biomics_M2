# -*- coding: utf-8 -*-
# arg1 : reference fasta file
# arg2 : bed file
# arg3 : int : window size (optionnal : default 101)

import matplotlib.pyplot as plt
import sys
from sequana import *


filename_ref = str(sys.argv[1])
filename_bed = str(sys.argv[2])

if len(sys.argv) > 3:
	window = int(sys.argv[3])
	# if even number, add 1 to get odd number
	if not window % 2:
		window += 1
else:
	window = 101


cov = GenomeCov(filename_bed)
# compute gc content with window
print("Compute GC content")
cov.compute_gc_content(filename_ref,circular = True, window_size=window)

all_opt_window = []

for chromosome in cov.chr_list:
	print("Compute optimal window for chromosome %s" % chromosome.chrom_name )
	#chromosome.df.sum()
	opt_window = int(round(chromosome.get_max_gc_correlation(filename_ref)))
	# if even number, add 1 to get odd number
	if not opt_window % 2:
		opt_window += 1
	all_opt_window.append(opt_window)
	plt.title("Optimal window = %d" % opt_window)
	plt.savefig('Cor_GC_coverage_window_size_%s.png' % chromosome.chrom_name)
	plt.show()
	plt.close()
	#chromosome.df.sum()

print(all_opt_window)
# compute gc content with optimal window

for opt in all_opt_window:
	cov.compute_gc_content(filename_ref,circular = True, window_size=opt)
	print("ok")
	for chromosome in cov.chr_list:
		chromosome.plot_gc_vs_coverage()
		plt.title("%s - Optimal window = %d" % (chromosome.chrom_name,opt_window))
		plt.savefig('Cor_GC_coverage_optimal_window_%s.png' % chromosome.chrom_name )
		plt.show()
		plt.close()



