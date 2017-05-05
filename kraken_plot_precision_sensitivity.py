# -*- coding: utf-8 -*-
# arg1 : fof_classification_result.txt : fof of classification results (csv files, ouput of kraken_genus_accuracy.py)
# arg2 : labels.txt : labels to use for each csv file, in the same order than fof_classification_result.txt
# arg3 : plot.png : name of plot. if "show" don't plot, show instead

# plot Sensitivity vs Precision of kraken results

# columns of input csv files must contain :
# "good_classification_at_level","wrong_classification_at_level","unknown_taxon_at_level",
# "good_classification_above_level","wrong_classification_above_level","unknown_taxon_above_level","Unclassified"

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
	xlim = [0.995,1.]
	ylim = [0.33,1.]

colors = ["b","b","b","b","b","b","k","k","k","k","r"]
marker = ["o","o","o","o","o","o","o","X",".","s","X"]


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



pylab.figure(figsize=(5, 5))

for i in range(len(list_input)):
	df = pd.read_csv(list_input[i])
	good_class_rank = df["good_classification_at_level"].sum()
	tot_class_rank = df["good_classification_at_level"].sum() + df["wrong_classification_at_level"].sum()
	tot = df["Unclassified"].sum() + tot_class_rank + df["unknown_taxon_at_level"].sum()
	#print(tot_class_rank)
	wrong_class_above = df["wrong_classification_above_level"].sum()
	#s = good_class_rank / float(tot_class_rank)
	s = good_class_rank / float(tot)
	p = good_class_rank / float(tot_class_rank + wrong_class_above)
	print("%s : Precision = %f, Sensitivity = %f" %(list_labels[i], p, s))
	pylab.plot(p,s,colors[i]+marker[i],label=list_labels[i])


pylab.title(title_plot)
pylab.legend(loc=2)
if manual_scale:
	pylab.xlim(xlim)
	pylab.ylim(ylim)
pylab.xlabel(x_label)
pylab.ylabel(y_label)
if save_plot:
	pylab.savefig(filename_plot)
else:
	pylab.show()







