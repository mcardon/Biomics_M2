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

M = 200 # param fft
alpha = 0.5
figsize = (8,10)
################################ EXECUTE ##############################################################################################

seq = DNA(filename_fasta)
seq.window = 0.05

len_genome = len(seq)

seq._myfft_gc_skew(M)

rep = Repeats(filename_fasta) 

################################ PLOT ##############################################################################################


x = range(len_genome)

fig, axarr = pylab.subplots(5,1, sharex=True, figsize=figsize)

############ GC skew
axarr[0].plot(x, (seq.AT_skew).T, 'b-', alpha=alpha)
# Make the y-axis label, ticks and tick labels match the line color.
axarr[0].set_ylabel('AT skew', color='k')

ax2 = axarr[0].twinx()
ax2.plot(x, ((seq.AT_skew).T).cumsum(), 'b--', alpha=alpha)
ax2.set_ylabel('AT skew\n(cumulative)', color='k')


############ AT skew
axarr[1].plot(x, (seq.GC_skew).T, 'b-', alpha=alpha)
axarr[1].set_ylabel('GC skew', color='k')

ax3 = axarr[1].twinx()
ax3.plot(x, ((seq.GC_skew).T).cumsum(), 'b--', alpha=alpha)
ax3.set_ylabel('GC skew\n(cumulative)', color='k')


############ Z curve
axarr[2].plot(x, seq._Xn/np.std(seq._Xn), 'k--', alpha=alpha)
axarr[2].plot(x, seq._Yn/np.std(seq._Yn), 'k-', alpha=alpha)
axarr[2].plot(x, seq._Zn/np.std(seq._Zn), 'k:', alpha=alpha)
axarr[2].set_ylabel('Z curve', color='k')

############ Z curve
axarr[3].plot(x, seq._c_fft, 'g-', alpha=alpha)
axarr[3].set_ylabel('FFT', color='k')

pylab.show()