# -*- coding: utf-8 -*-
# arg1 : blasr output file (with output option -m 5 : should be 19 columns)
# arg2 : optionnal : save, saveonly : (save : will save table of results) (saveonly : save table, don't plot anything)


################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from random import shuffle
import sys

################################ PARAMETERS ##############################################################################################

file_blasr = str(sys.argv[1])

save_result = False
do_plots = True
if len(sys.argv) > 2:
	if str(sys.argv[2]) == "save":
		save_result = True
	elif str(sys.argv[2]) == "saveonly":
		save_result = True
		do_plots = False

# file format 
# https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format
blasr_columns = "qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq".split(" ")
colormap = "gist_rainbow"

threshold_variants = 50

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
	# extract alignment columns
	res_align = res_blasr[["qName","qLength","qAlignedSeq","tAlignedSeq","matchPattern"]]
	# drop alignment columns (useless and too heavy)
	# res_blasr.drop("qAlignedSeq",axis=1)
	res_blasr.drop(["qAlignedSeq","tAlignedSeq","matchPattern"], axis=1, inplace=True)

	return res_blasr, res_align


################################ IMPORT DATA ##############################################################################################

# read file 
res_blasr, res_align = read_blasr_result(file_blasr, blasr_columns)


################################ EXECUTE ##############################################################################################

# normalise score by length
res_blasr["score_norm"] = res_blasr["score"] / res_blasr["qLength"]

nb_contigs = res_align.shape[0]
df_variants = pd.DataFrame()
# for each contig, look the size of variations (rough caracterisation : will be better with bam output from blasr when available)
ll_size_variants = []
ll_insertions = []
ll_deletions = []
for ind_contig in range(nb_contigs):
	name = res_align.loc[ind_contig,"qName"]
	match_pattern = res_align.loc[ind_contig,"matchPattern"]
	seq_q = res_align.loc[ind_contig,"qAlignedSeq"]
	seq_t = res_align.loc[ind_contig,"tAlignedSeq"]

	seq_q = "".join([seq_q[i] if match_pattern[i] == "*" else "|" for i in range(len(seq_q)) ])
	seq_t = "".join([seq_t[i] if match_pattern[i] == "*" else "|" for i in range(len(seq_t)) ])

	### save length of variants
	# get rid of matches
	insertions = seq_q.split("|")
	deletions = seq_t.split("|")

	# filter out empty strings
	insertions = [ variant_i for variant_i in insertions if variant_i]
	deletions = [ variant_i for variant_i in deletions if variant_i]

	# size of variants
	variants_size = [ len(variant_i) for variant_i in insertions if variant_i]

	# filter out mismatches
	not_mismatch = [i for i in range(len(insertions)) if not (bool(insertions[i].replace('-','')) & bool(deletions[i].replace('-','')))]
	insertions = [insertions[i] for i in not_mismatch if len(insertions[i].replace('-','')) > threshold_variants ]
	deletions = [deletions[i] for i in not_mismatch if len(deletions[i].replace('-','')) > threshold_variants ]

	# append result
	ll_insertions.append(insertions)
	ll_deletions.append(deletions)
	ll_size_variants.append(variants_size)

	# convert results to df
	df_to_append = pd.DataFrame()
	df_to_append["seq_indel"] = insertions + deletions
	df_to_append["type_indel"] = ["insertion"]*len(insertions) + ["deletion"]*len(deletions)
	df_to_append["contig"] = [name]*df_to_append.shape[0]

	# append result df
	df_variants = df_variants.append(df_to_append)

# set index
df_variants.index = range(df_variants.shape[0])




################################ SAVE ##############################################################################################

# save clean result
if save_result:
	res_blasr.to_csv(file_blasr + "_scores.csv",index=False)
	df_variants.to_csv(file_blasr + "_variants.csv",index=False)


################################ PLOTS ##############################################################################################


if do_plots:
	# title
	title_plot = file_blasr.split("/")[-1]
	# colors
	cmap = pylab.cm.get_cmap(colormap)
	# shuffel colors :  in case 2 adjacent contigs have the same color, user can plot again to see better
	shuffle_col = list(np.linspace(0,1,res_blasr.shape[0]))
	shuffle(shuffle_col)
	colors = [cmap(i) for i in shuffle_col]

	# plot score nor normalised
	pylab.plot(res_blasr["qLength"], res_blasr["score"], "bo" ,alpha=0.5)
	pylab.xlabel("Length of contig")
	pylab.ylabel("Score blasr (not normalised)")
	pylab.title(title_plot)
	pylab.show()

	# plot score normalised by lenght
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


	# ###### histograms length of variants
	# fig, axarr = pylab.subplots(nb_contigs,figsize=(12,2*nb_contigs), sharex=True)
	# #fig.suptitle("Size of variants by contig", fontsize=14)

	# bin_to_use = max([max(variants_size) if len(variants_size) > 0 else 0 for variants_size in ll_size_variants ])

	# for ind_contig in range(nb_contigs):
	# 	ax = axarr[ind_contig]
	# 	# plot coverage found by blasr, with score
	# 	variants_size = ll_size_variants[ind_contig]
	# 	insertions = ll_insertions[ind_contig]
	# 	deletions = ll_deletions[ind_contig]
	# 	ax.hist(variants_size, label=str(ind_contig), bins=int(bin_to_use/10))
	# 	ax.set_title(res_align.loc[ind_contig,"qName"] + ' , Length : ' + str(res_align.loc[ind_contig,"qLength"]) + ', total variants : ' + str(sum(variants_size)))
	# 	ax.set_xlim([-10, bin_to_use+10])
	# 	ax.set_yscale('log')
	# 	ax.set_ylabel("log(#)")

	# ax.set_xlabel("Size of variants")
	# fig.tight_layout()
	# #pylab.savefig("./../test_fig_XXX.png")
	# pylab.show()
	# pylab.close("all")














