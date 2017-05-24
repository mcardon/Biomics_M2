# -*- coding: utf-8 -*-
# arg1 : file of filenames : list of bam file (absolute path)
# arg2 : csv file description
# arg3 : analysis name (for output files and plots) : if show, don't save, show instead

################################ IMPORT ################################################################################################
import os
import sys

from sequana import pacbio

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

figsize_read_len = (6,4)
figsize_GC       = (6,4)
figsize_ZMW      = (6,4)
figsize_SNR      = (6,4)

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
	list_BAM.append(pacbio.BAMPacbio(f))



# plot read length
fig, ax = pylab.subplots(1,1, figsize=figsize_read_len)
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	bam.hist_len(hold=True,grid=False,label=labels[i],title="")
ax.legend()
fig.tight_layout()
if save_plots:
	pylab.savefig(filename_output.replace(".","_read_len."), dpi=182)
	pylab.clf()
else:
	pylab.show()


# plot GC %
fig, ax = pylab.subplots(1,1, figsize=figsize_GC)
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	bam.hist_GC(hold=True,label=labels[i],grid=False,title="")
# pylab.title("GC content")
pylab.legend()
fig.tight_layout()
if save_plots:
	pylab.savefig(filename_output.replace(".","_GC."), dpi=182)
	pylab.clf()
else:
	pylab.show()

# plot ZMW passes
fig, ax = pylab.subplots(1,1, figsize=figsize_ZMW)
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	bam.hist_ZMW_subreads(hold=True,label=labels[i],title="",grid=False,xlabel="Number of passes")
	pylab.xlim([0,45])

	pylab.legend()
	fig.tight_layout()
	if save_plots:
		pylab.savefig(filename_output.replace(".","_ZMW_%d."%i), dpi=182)
		pylab.clf()
	else:
		pylab.show()


# plot snr
for i in range(len(list_BAM)):
	bam = list_BAM[i]
	# plot read length
	fig, ax = pylab.subplots(1,1, figsize=figsize_SNR)
	bam.hist_snr(grid=False,title="")
	pylab.xlim([3,16])
	# pylab.title("SNR %s" %labels[i])
	pylab.legend()
	fig.tight_layout()
	if save_plots:
		pylab.savefig(filename_output.replace(".","_SNR_%d." %i), dpi=182)
		pylab.clf()
	else:
		pylab.show()



