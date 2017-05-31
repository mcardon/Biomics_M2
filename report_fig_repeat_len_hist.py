# -*- coding: utf-8 -*-
# arg1 : reference_genome.fasta
# arg2 : filename for output plot (if "show" dont save, show instead)

################################ IMPORT ################################################################################################
import pylab
from sequana import sequence
import collections
import sys

################################ PARAMETERS ##############################################################################################

f_input         = str(sys.argv[1])
filename_output = str(sys.argv[2])

save_output = True
if filename_output == "show":
	save_output = False


alpha = 0.6
nb_bars = 30
figsize = (7,7)
max_chars = 100
fontsize = 16

if save_output:
	title = filename_output.replace(".png","").replace("_"," ").replace(".txt","")
else:
	title = f_input.replace("fof_","").replace(".txt","")


################################ INPUT DATA ##############################################################################################
rep = sequence.Repeats(f_input)
rep.threshold = 100



fig, ax = pylab.subplots(1,1, figsize=figsize)
rep.hist_length_repeats(bins=nb_bars, alpha=alpha, fontsize=fontsize,
 grid=False, title="Repeat length", xlabel="Repeat length", ylabel="#")
pylab.tight_layout()

if save_output:
	pylab.savefig(filename_output,dpi=96)
	pylab.clf()
else:
	pylab.show()


"""
nb_bars = 30
rep.hist_length_repeats(bins=nb_bars, alpha=alpha, fontsize=fontsize,
 grid=False, title="Repeat length", xlabel="Repeat length", ylabel="#")
pylab.show()
"""
