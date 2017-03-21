# -*- coding: utf-8 -*-
# arg1 : blasr output file (with output option -m 5 : should be 19 columns)
# arg2 : optionnal : if save : will save table of results


################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np

import sys

################################ PARAMETERS ##############################################################################################

file_blasr = str(sys.argv[1])

save_result = False
if len(sys.argv) > 2:
	if str(sys.argv[1]) == "save":
		save_result = True

# file format 
# https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
blasr_columns = "qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq".split(" ")

################################ FUNCTIONS ##############################################################################################

def read_blasr_result(file_blasr, blasr_columns):
	res_blasr = []
	f = open(file_blasr, 'r')
	l = f.readline()
	while(l):
		contig = l.replace("\n","").replace("  "," ").split(" ")
		res_blasr.append(contig)
		l = f.readline()
	f.close()

	res_blasr = pd.DataFrame(res_blasr)
	res_blasr.columns = blasr_columns

	# convert to numeric
	res_blasr["qLength"] = res_blasr["qLength"].astype(int)
	res_blasr["score"] = res_blasr["score"].astype(float)
	# drop alignment columns (useless and too heavy)
	res_blasr.drop("qAlignedSeq",axis=1)
	res_blasr.drop(["qAlignedSeq","tAlignedSeq","matchPattern"], axis=1, inplace=True)

	return res_blasr


################################ IMPORT DATA ##############################################################################################

# read file 
res_blasr = read_blasr_result(file_blasr, blasr_columns)

# normalise score by length
res_blasr["score_norm"] = res_blasr["score"] / res_blasr["qLength"]

# save clean result
if save_result:
	res_blasr.to_csv(file_blasr + "_scores.csv")



################################ PLOTS ##############################################################################################

title_plot = file_blasr.split("/")[-1]

pylab.plot(res_blasr["qLength"], res_blasr["score"],"bo",alpha=0.5)
pylab.xlabel("Length of contig")
pylab.ylabel("Score blasr (not normalised)")
pylab.title(title_plot)
pylab.show()


pylab.plot(res_blasr["qLength"], res_blasr["score_norm"],"bo",alpha=0.5)
pylab.xlabel("Length of contig")
pylab.ylabel("Score blasr (normalised by length)")
pylab.title(title_plot)
pylab.show()

# plot coverage found by blasr, with score
pylab.figure(figsize=(15,3))
for i in range(res_blasr.shape[0]):
	contig = res_blasr["qName"][i]
	res_to_plot = res_blasr[res_blasr["qName"] == contig]
	start = int(res_to_plot["tStart"].values[0])
	end = int(res_to_plot["tEnd"].values[0])
	score = float(res_to_plot["score_norm"].values[0])
	pylab.plot([start, end],[score]*2, "b-" )
pylab.title("Coverage by contigs (blasr)\n%s" % title_plot)
pylab.xlabel("Reference genome position")
pylab.ylabel("Score blasr (normalised by length)")
pylab.tight_layout()
pylab.show()




