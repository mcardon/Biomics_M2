# -*- coding: utf-8 -*-
# arg1 : blasr output file (with output option -m 5 : should be 19 columns)
# arg2 : optionnal : if save : will save table of results


################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from random import shuffle
import sys

################################ PARAMETERS ##############################################################################################

file_blasr = str(sys.argv[1])

save_result = False
if len(sys.argv) > 2:
	if str(sys.argv[2]) == "save":
		save_result = True

# file format 
# https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
blasr_columns = "qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq".split(" ")
colormap = "gist_rainbow"

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

# title
title_plot = file_blasr.split("/")[-1]
# colors
cmap = pylab.cm.get_cmap(colormap)
# shuffel colors :  in case 2 adjacent contigs have the same color, user can plot again to see better
shuffle_col = list(np.linspace(0,1,res_blasr.shape[0]))
shuffle(shuffle_col)
colors = [cmap(i) for i in shuffle_col]



pylab.plot(res_blasr["qLength"], res_blasr["score"], "bo" ,alpha=0.5)
pylab.xlabel("Length of contig")
pylab.ylabel("Score blasr (not normalised)")
pylab.title(title_plot)
pylab.show()


pylab.plot(res_blasr["qLength"], res_blasr["score_norm"],"bo",alpha=0.5)
pylab.xlabel("Length of contig")
pylab.ylabel("Score blasr (normalised by length)")
pylab.title(title_plot)
pylab.show()


fig, axarr = pylab.subplots(2,figsize=(15,8), sharex=True)
fig.suptitle("Coverage by contigs (blasr)\n%s" % title_plot, fontsize=10)
# plot coverage found by blasr, with score
ax = axarr[0]
for i in range(res_blasr.shape[0]):
	res_to_plot = res_blasr.iloc[i,:]
	contig = res_to_plot["qName"]
	start = int(res_to_plot["tStart"])
	end = int(res_to_plot["tEnd"])
	score = float(res_to_plot["score_norm"])
	ax.plot([start, end],[score]*2, ls='-', lw=5, color=colors[i], solid_capstyle="butt" )
ax.set_ylabel("Score blasr (normalised by length)")



# plot coverage found by blasr, with random y distribution (to see if there are overlaps)
ax = axarr[1]
y = list(np.linspace(0,1,res_blasr.shape[0]))
for i in range(res_blasr.shape[0]):
	res_to_plot = res_blasr.iloc[i,:]
	contig = res_to_plot["qName"]
	start = int(res_to_plot["tStart"])
	end = int(res_to_plot["tEnd"])
	score = float(res_to_plot["score_norm"])
	ax.plot([start, end],[y[i]]*2, ls='-', lw=5, color=colors[i], solid_capstyle="butt" )
ax.set_ylabel("Random")
ax.set_xlabel("Reference genome position")


fig.subplots_adjust(bottom=0.2,top=0.8)
#fig.tight_layout()
pylab.show()
pylab.close("all")
