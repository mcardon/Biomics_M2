# -*- coding: utf-8 -*-
# arg1 : file of filenames : list of bam file (absolute path)
# arg2 : csv file description
# arg3 : analysis name (for output files and plots) : if show, don't save, show instead

################################ IMPORT ################################################################################################
import os
import sys

scriptpath = "/home/mcardon/Mel/Code/sequana/sequana/"
sys.path.append(os.path.abspath(scriptpath))
import pacbio

from sequana.lazy import pylab
from sequana.lazy import pandas as pd

################################ PARAMETERS ##############################################################################################

fof_BAM = str(sys.argv[1])
file_description = str(sys.argv[2])
filename_output = str(sys.argv[3])

if filename_output == "show":
	save_plots = False
else:
	save_plots = True

################################ FUNCTIONS ##############################################################################################

def read_fof(fof_BAM):
	# list of all files
	f = open(fof_BAM, 'r')
	list_files = []
	for name in f.readlines():
		f_name = name.split('\n')[0]
		list_files.append(f_name)
	f.close()
	return list_files


################################ IMPORT DATA ##############################################################################################

list_files = read_fof(fof_BAM)
description = pd.read_csv(file_description)

list_BAM = []
labels = []
for f in list_files:
	short_name_bam = f.split("/")[-1]
	labels.append(description[description["Filename"] == short_name_bam]["polymerase"].values[0])
	list_BAM.append(pacbio.PacBioInputBAM(f))



# plot read length
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	bam.hist_len(hold=True,label=labels[i])
pylab.title("Read length")
pylab.legend()
pylab.show()

# plot GC %
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	# plot read length
	bam.hist_GC(hold=True,label=labels[i])
pylab.title("GC content")
pylab.legend()
pylab.show()

# plot ZMW passes
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	# plot read length
	bam.hist_ZMW_subreads(hold=True,label=labels[i])
pylab.title("ZMW passes")
pylab.legend()
pylab.show()

# plot snr
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	# plot read length
	bam.hist_snr()
	pylab.title("SNR %s" %labels[i])
	pylab.legend()
	pylab.show()



