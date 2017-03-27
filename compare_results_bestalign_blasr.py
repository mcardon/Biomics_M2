# -*- coding: utf-8 -*-
# arg1 : fof.txt : file of filenames : csv files with normalised blasr score (as returned by compare_contigs_blasr.py)
# corresponds to multiple alignements of one assembly : evaluate quality of assembly with multiple references
# arg2 : output.csv (if no, will not save output)
# arg3 : optionnal : show or filename.png : if show : show figure, else save figure as filename.png (if nothing, no plot)
# arg4 : optionnal : filename_not_covered.csv : if any, save the list of not covered regions : must show or save plot to do this

## file format : /home/mcardon/Mel/Datas/2017_03_17_De_novo_smrtlink/blasr_results_multiple_ref/100X_ref_NC_rotate_1225857.fasta.blasr5.out_scores.csv
## string between 002929_ and .fasta will be parsed : must contain "rotate" and a number for the rotation step


################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from random import shuffle
import sys


################################ PARAMETERS ##############################################################################################

fof_scores_csv = str(sys.argv[1])
file_output = str(sys.argv[2])

save_csv = True
if file_output == 'no':
	save_csv = False

do_plot = False
if len(sys.argv) > 3:
	file_plot = str(sys.argv[3])
	do_plot = True
	title_plot = "25X De novo Pacbio"
	colormap = "gist_rainbow"

	if file_plot == 'show':
		save_plot = False
	else:
		save_plot = True

save_not_covered = False
if len(sys.argv) > 4:
	filename_not_covered = str(sys.argv[4])
	save_not_covered = True

len_genome = 4086189

################################ IMPORT DATA ##############################################################################################

# list of all files
f = open(fof_scores_csv, 'r')
list_files = []
for name in f.readlines():
	f_name = name.split('\n')[0]
	list_files.append(f_name)
f.close()

# concat all result files
list_to_concat = []
for file_result in list_files:
	df = pd.read_csv(file_result,sep=",")
	short_name = file_result.split("/")[-1].split("002929_")[-1].split(".fasta")[0]
	short_name = short_name.split("_")
	df["ref_type"] = [short_name[0]]*df.shape[0]
	if len(short_name) > 1:
		df["rotation"] = [int(short_name[1])]*df.shape[0]
	else:
		df["rotation"] = [0]*df.shape[0]
	list_to_concat.append(df)
#list_to_concat = [pd.read_csv(file_result,sep=",") for file_result in list_files]
result = pd.concat(list_to_concat)


################################ EXECUTE ##############################################################################################

# for each contig, choose the alignement with best score
set_contigs = set(result["qName"])

best_contigs = []
for contig in set_contigs:
	df = result[result["qName"] == contig]
	best_score = min(df["score_norm"])
	#print(contig, best_score)
	best_contigs.append(pd.DataFrame(df[df["score_norm"] == best_score].iloc[0,:]).T)

res_best = pd.concat(best_contigs)
if not isinstance(res_best, pd.DataFrame):
	res_best = pd.DataFrame(res_best).T

if save_csv:
	res_best.to_csv(file_output,index=False)

print("Best score is %.2f" %(round(min(res_best["score_norm"]),2)))
print("Worst score is %.2f" %(round(max(res_best["score_norm"]),2)))


################################ FUNCTIONS FOR PLOT ##################################################################################

def plot_contigs(res_best, ax, mode="score"):
	"""
	Plot contig aligned on genome with given y position
	mode : string : score or random : choose y position
	"""

	# create y array if needed (plot without score)
	if mode != "score":
		y = list(np.linspace(0,1,res_best.shape[0]))

	# save result of contigs
	list_contigs = []

	# plot contigs
	for i in range(res_best.shape[0]):
		res_to_plot = res_best.iloc[i,:]
		contig = res_to_plot["qName"]
		start = (int(res_to_plot["tStart"]) + int(res_to_plot["rotation"]))
		end = (int(res_to_plot["tEnd"]) + int(res_to_plot["rotation"]))
		# check if start and end need to be rotated
		if int(abs(start/len_genome)) == 1:
			start = start % len_genome
			end = end % len_genome
		# save contigs
		list_contigs.append([start,end])
		# choose y postion
		if mode == "score":
			score = float(res_to_plot["score_norm"])
			ax.plot([start, end],[score]*2, ls='-', lw=5, color=colors[i], solid_capstyle="butt" )
		else:
			ax.plot([start, end],[y[i]]*2, ls='-', lw=5, color=colors[i], solid_capstyle="butt" )

	return list_contigs


def areas_not_covered(list_contigs, len_genome):
	"""
	return list of areas not covered by contigs
	"""
	# plot areas not covered
	start_first_contig = min([contig[0] for contig in list_contigs])
	#end_last_contig = max([contig[1] for contig in list_contigs])
	genome_not_covered = [[start_first_contig, start_first_contig + len_genome]]
	for contig in list_contigs:
		start = contig[0]
		end = contig[1]

		prev_genome_part = [0,0]
		i = 0
		inside_contig = False
		while(i < len(genome_not_covered)):
			genome_part = genome_not_covered[i]

			if inside_contig:
				if (end < genome_part[0]):
					# end before next uncovered region
					genome_not_covered = new_genome_not_covered + genome_not_covered[i:len(genome_not_covered)]
					inside_contig = False
					#print("should break here 3")
					break
				elif (end >= genome_part[0]) & (end <= genome_part[1]):
					genome_not_covered = new_genome_not_covered + [[end+1, genome_part[1]]] + genome_not_covered[min(i+1, len(genome_not_covered)):len(genome_not_covered)]
					inside_contig = False
					#print("should break here 2")
					break
			else:
				# whole contig in already covered area -> do nothing
				if (start > prev_genome_part[1]) & (end < genome_part[0]):
					#print("should break here 4")
					break
				# start inside a non covered area
				if (start >= genome_part[0]) & (start <= genome_part[1]):
					inside_contig = True
					# remove this part of genome, replace by not covered part
					new_genome_not_covered = genome_not_covered[0:i] + [[genome_part[0], max(start-1, genome_part[0])]]
					#print("start found")
					# is the end inside this part ?
					if (end <= genome_part[1]):
						genome_not_covered = new_genome_not_covered + [[end+1, genome_part[1]]] + genome_not_covered[min(i+1, len(genome_not_covered)):len(genome_not_covered)]
						#print("should break here 1")
						inside_contig = False
						break
			prev_genome_part = genome_not_covered[i]
			#print(i, genome_not_covered[i])
			i +=1

		# if still inside contig, need to finish outside loop
		if inside_contig:
			# print("here")
			# print(contig)
			# print(new_genome_not_covered)
			# print(genome_not_covered)
			genome_not_covered = new_genome_not_covered

	# clean result
	genome_not_covered = [area for area in genome_not_covered if area[0] < area[1]]

	return genome_not_covered


################################ PLOT ##############################################################################################

if do_plot:

	cmap = pylab.cm.get_cmap(colormap)
	# shuffle colors :  in case 2 adjacent contigs have the same color, user can plot again to see better
	shuffle_col = list(np.linspace(0,1,res_best.shape[0]))
	shuffle(shuffle_col)
	colors = [cmap(i) for i in shuffle_col]
	
	pylab.plot(res_best["qLength"], res_best["score_norm"],"bo",alpha=0.5)
	pylab.xlabel("Length of contig")
	pylab.ylabel("Score blasr (normalised by length)")
	pylab.title(title_plot)
	if save_plot:
		pylab.savefig(file_plot.replace(".png","_scores.png"))
	else:
		pylab.show()

	fig, axarr = pylab.subplots(2,figsize=(15,8), sharex=True)
	fig.suptitle("Coverage by contigs (blasr)\n%s" % title_plot, fontsize=10)
	# plot coverage found by blasr, with score
	ax = axarr[0]
	list_contigs = plot_contigs(res_best, ax, mode="score")
	genome_not_covered = areas_not_covered(list_contigs, len_genome)
	# add grey on not covered areas
	for area in genome_not_covered:
		ax.axvspan(area[0],area[1], alpha=0.1, color='k')
	ax.set_ylabel("Score blasr (normalised by length)")

	# plot coverage found by blasr, with random y distribution (to see if there are overlaps)
	ax = axarr[1]
	list_contigs = plot_contigs(res_best, ax, mode="random")
	for area in genome_not_covered:
		ax.axvspan(area[0],area[1], alpha=0.1, color='k')

	ax.set_ylabel("Random")
	ax.set_xlabel("Reference genome position")

	fig.subplots_adjust(bottom=0.2,top=0.8)
	#fig.tight_layout()
	if save_plot:
		pylab.savefig(file_plot)
	else:
		pylab.show()
	pylab.close("all")

	if save_not_covered:
		df = pd.DataFrame(genome_not_covered)
		df.to_csv(filename_not_covered, index=None,header=False)





