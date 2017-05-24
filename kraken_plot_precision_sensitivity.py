# -*- coding: utf-8 -*-
# arg1 : fof_classification_result.txt : fof of classification results (csv files, ouput of kraken_genus_accuracy.py)
# arg2 : labels.txt : labels to use for each csv file, in the same order than fof_classification_result.txt
# arg3 : plot.png : name of plot. if "show" don't plot, show instead

# plot Sensitivity vs Precision of kraken results

# columns of input csv files must contain :
# "good_classification_at_level","wrong_classification_at_level","unknown_taxon_at_level",
# "good_classification_above_level","wrong_classification_above_level","unknown_taxon_above_level","Unclassified","total_N_reads"

# edit parameters for custom plots

################################ IMPORT ################################################################################################
#import os
import sys
from sequana.lazy import pylab
from sequana.lazy import pandas as pd

################################ PARAMETERS ##############################################################################################

fof_input     = str(sys.argv[1])
labels_input  = str(sys.argv[2])
filename_plot = str(sys.argv[3])

save_plot = True
if filename_plot == "show":
	save_plot = False

level_tax = "Genus"
title_plot = "Classification at %s level" % level_tax
x_label = "%s Precision" %level_tax
y_label = "%s Sensitivity"%level_tax

manual_scale = True
if manual_scale:
	xlim = [0.86,1.]
	ylim = [0.,1.]

# colors = ["b","k","r","g","y","m","c"]
# colors = ["k","r","g","y","m","c"]
colors = ["r","g","y","m","c"]
marker = ["o","o","o","o","o","o","o"]

loc_label = 3
################################ INPUT DATA ##############################################################################################

# list of input files
f = open(fof_input, 'r')
list_input = [name.split('\n')[0] for name in f.readlines()]
f.close()

# list of labels
f = open(labels_input, 'r')
list_labels = [name.split('\n')[0] for name in f.readlines()]
f.close()


################################ EXECUTE ##############################################################################################




#pylab.figure(figsize=(5, 5))
fig1, ax1 = pylab.subplots(1,1, figsize=(5, 5))
fig2, ax2 = pylab.subplots(1,1, figsize=(5, 5))

res_PS = []

for i in range(len(list_input)):
	df = pd.read_csv(list_input[i])

	##### Precision without unknown taxons
	good_class_rank = df["good_classification_at_level"].sum()
	tot_class_rank = df["good_classification_at_level"].sum() + df["wrong_classification_at_level"].sum()
	tot = df["total_N_reads"].sum()
	wrong_class_above = df["wrong_classification_above_level"].sum()

	##### Precision without unknown taxons
	# some reads are classified but we dont find any info : cannot ignore them
	# wrong_class_above_u = df["wrong_classification_above_level"].sum() + df["unknown_taxon_above_level"].sum()
	# wrong_not_output = tot - (tot_class_rank + df["good_classification_above_level"].sum() + wrong_class_above_u )
	# wrong_class_above_u += wrong_not_output

	wrong_class_above_u = tot - (good_class_rank + df["good_classification_above_level"].sum() + df["Unclassified"].sum())
	s = good_class_rank / float(tot)
	p = good_class_rank / float(tot_class_rank + wrong_class_above)
	p_u = good_class_rank / float(tot_class_rank + wrong_class_above_u)
	print("%s : Precision = %f, Precision with unknowns = %f, Sensitivity = %f" %(list_labels[i], p, p_u, s))
	ax1.plot(p,s,colors[i]+marker[i],label=list_labels[i])
	ax2.plot(p_u,s,colors[i]+marker[i],label=list_labels[i])
	res_PS.append([list_labels[i], s, p, p_u])


fig1.suptitle(title_plot)
fig2.suptitle(title_plot)
ax1.legend(loc=loc_label)
ax2.legend(loc=loc_label)

if manual_scale:
	ax1.set_xlim(xlim)
	ax1.set_ylim(ylim)
	ax2.set_xlim(xlim)
	ax2.set_ylim(ylim)

ax1.set_xlabel(x_label)
ax1.set_ylabel(y_label)
ax2.set_xlabel(x_label)
ax2.set_ylabel(y_label)

fig1.tight_layout()
fig1.subplots_adjust(top=0.92)
fig2.tight_layout()
fig2.subplots_adjust(top=0.92)

if save_plot:
	fig1.savefig(filename_plot.replace(".png","") + "_without_unknown.png")
	fig2.savefig(filename_plot.replace(".png","") + "_with_unknown.png")
	df_res_PS = pd.DataFrame(res_PS)
	df_res_PS.columns = ["label","sensitivity","precision","precison_with_unknown"]
	df_res_PS.to_csv(filename_plot.replace(".png","") + "_data.csv",index=None)
else:
	fig1.show()
	fig2.show()





