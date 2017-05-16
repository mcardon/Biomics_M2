# -*- coding: utf-8 -*-
# arg1 : filename genome (fasta format)
# arg2 : prefix filename for plots (if show : don't save anything, show plots)

################################ IMPORT ################################################################################################
import pylab
import numpy as np
from sequana.sequence import DNA, Repeats

import sys

################################ PARAMETERS ##############################################################################################

filename_fasta = str(sys.argv[1])
filename_output = str(sys.argv[2])

save_output = True
if filename_output == "show":
	save_output = False

# sliding window for GC skew
window = 0.1

# threshold to use to output repeats
threshold = 100

# what to do
do_skew = True
do_CDS = False
do_rep_histogram = False
do_FFT = False

M = 2000 # param fft
alpha = 0.5
figsize = (8,10)
fontsize = 16
title = "Bordetella Bronseptica genome"

################################ EXECUTE ##############################################################################################

seq = DNA(filename_fasta)

len_genome = len(seq)

if do_rep_histogram:
	rep = Repeats(filename_fasta)
	rep.threshold = threshold
	rep.do_merge = False

if do_skew:
	seq.window = window

	if do_FFT:
		seq._myfft_gc_skew(M)

	rep = Repeats(filename_fasta)
	rep.threshold = threshold
	rep.do_merge  = True


################################ PLOT ##############################################################################################

if do_skew:
	x = range(len_genome)
	if do_FFT:
		fig, axarr = pylab.subplots(4,1, sharex=True, figsize=figsize, gridspec_kw = {'height_ratios':[3,3,3,1]})
	else:
		fig, axarr = pylab.subplots(3,1, sharex=True, figsize=figsize)

	############ GC skew
	axarr[0].plot(x, (seq.AT_skew).T, 'b-', alpha=alpha)
	for repeat in rep.begin_end_repeat_position:
		axarr[0].axvspan(repeat[0], repeat[1], alpha=0.1, color='black')
	# Make the y-axis label, ticks and tick labels match the line color.
	axarr[0].set_ylabel('AT skew', color='k', fontsize=fontsize)

	ax2 = axarr[0].twinx()
	ax2.plot(x, ((seq.AT_skew).T).cumsum(), 'b--', alpha=alpha)
	for repeat in rep.begin_end_repeat_position:
		axarr[1].axvspan(repeat[0], repeat[1], alpha=0.1, color='black')
	ax2.set_ylabel('AT skew\n(cumulative)', color='k',fontsize=fontsize)

	############ AT skew
	axarr[1].plot(x, (seq.GC_skew).T, 'b-', alpha=alpha)
	axarr[1].set_ylabel('GC skew', color='k',fontsize=fontsize)

	ax3 = axarr[1].twinx()
	ax3.plot(x, ((seq.GC_skew).T).cumsum(), 'b--', alpha=alpha)
	ax3.set_ylabel('GC skew\n(cumulative)', color='k',fontsize=fontsize)

	############ Z curve
	axarr[2].plot(x, seq._Xn/np.std(seq._Xn), 'k--', alpha=alpha)
	axarr[2].plot(x, seq._Yn/np.std(seq._Yn), 'k-')
	axarr[2].plot(x, seq._Zn/np.std(seq._Zn), 'k:',alpha=alpha)
	axarr[2].set_ylabel('Z curve', color='k',fontsize=fontsize)

	if do_FFT:
		############ FFT
		axarr[3].plot(x, seq._c_fft, 'g-', alpha=alpha)
		axarr[3].set_ylabel('FFT', color='k',fontsize=fontsize)
		axarr[3].set_xlabel('Genome postion', color='k',fontsize=fontsize)
	else:
		axarr[2].set_xlabel('Genome postion', color='k',fontsize=fontsize)

	############ Repeats
	# for repeat in rep.begin_end_repeat_position:
	# 	axarr[4].axvspan(repeat[0], repeat[1], alpha=0.3, color='black')
	fig.suptitle(title, fontsize=(fontsize+2))

	pylab.tight_layout()
	fig.subplots_adjust(top=0.92)

	if save_output:
		pylab.savefig(filename_output, dpi=182)
		pylab.clf()
	else:
		pylab.show()

if do_CDS:
	pylab.clf()
	seq.barplot_count_ORF_CDS_by_frame()
	pylab.title("Bordetella Pertussis genome\nCount ORF and CDS by frame")
	if save_output:
		pylab.savefig(filename_output.replace(".","_barplot_CDS."), dpi=144)
		pylab.clf()
	else:
		pylab.show()

	seq.hist_ORF_CDS_logscale()
	pylab.title("Length distribution of ORF and CDS")
	if save_output:
		pylab.savefig(filename_output.replace(".","_hist_CDS_log."), dpi=144)
		pylab.clf()
	else:
		pylab.show()

if do_rep_histogram:
	pylab.clf()
	rep.hist_length_repeats(bins=110,grid=False,title="")
	if save_output:
		pylab.savefig(filename_output.replace(".","_hist_repeats."), dpi=144)
		pylab.clf()
	else:
		pylab.show()
