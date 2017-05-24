# -*- coding: utf-8 -*-
# arg1 : fof of files with data (Experiment_data.csv, output of kraken_plot_sensitivity.py)
# arg2 : prefix filename for plots (if show : don't save anything, show plots)

################################ IMPORT ################################################################################################
import pylab
from sequana.lazy import pandas as pd

import sys

################################ PARAMETERS ##############################################################################################

fof_input       = str(sys.argv[1])
filename_output = str(sys.argv[2])

save_output = True
if filename_output == "show":
	save_output = False

colors = ["y","c","b","g","r","m","c"]
# colors = ["k","r","g","y","m","c"]

## marker for line style
marker = ["--","--","--","--","--","o","o"]

marker_kmer_size = False
# if False, specify markers here for kmer size
marker_kmer = ["o","X","s","D","v","*","p"]

loc_label = 3
alpha = 0.6
manual_scale = True
if manual_scale:
	xlim = [0.84,1.]
	ylim = [0.45,1.]


min_kmer_size = 16
level_tax = "Genus"
x_label = "%s Precision" %level_tax
y_label = "%s Sensitivity"%level_tax

if save_output:
	title = filename_output.replace(".png","").replace("_"," ").replace(".txt","").replace("fof","")
else:
	title = fof_input.replace("fof_","").replace(".txt","")


################################ INPUT DATA ##############################################################################################

# list of input files
f = open(fof_input, 'r')
list_input = [name.split('\n')[0] for name in f.readlines()]
f.close()

df = []
for file_input in list_input:
	df.append(pd.read_csv(file_input))

df = pd.concat(df)
df["kmer_size"] = df["label"].str.split("_").str[-2]
df["kmer_size"] = df["kmer_size"].str.replace("k","").astype(int)

df["DB_size"] = df["label"].str.split("_").str[-1]

## filter out small kmer
df = df[df["kmer_size"] > min_kmer_size]

## plot lines
fig1, ax1 = pylab.subplots(1,1, figsize=(5, 5))
all_size = df["DB_size"].unique()
for i in range(len(all_size)):
	size = all_size[i]
	df_to_plot = df[df["DB_size"] == size]
	p = df_to_plot["precison_with_unknown"]
	s = df_to_plot["sensitivity"]
	ax1.plot(p,s,colors[i]+marker[i],label=size,alpha=alpha)

fig1.suptitle(title)


if marker_kmer_size:
	# plot legend before
	ax1.legend(loc=loc_label)

	## plot markers as kmer length
	all_kmer = df["kmer_size"].unique()
	for kmer in all_kmer:
		df_to_plot = df[df["kmer_size"] == kmer]
		p = df_to_plot["precison_with_unknown"]
		s = df_to_plot["sensitivity"]
		ax1.scatter(p, s, s=100, c="k", marker=r"$ {} $".format(kmer), edgecolors='none' )
else:
	## plot regular markers (see parameters)
	all_kmer = list(df["kmer_size"].unique())
	for i in range(len(all_kmer)):
		kmer = 	all_kmer[i]
		df_to_plot = df[df["kmer_size"] == kmer]
		p = df_to_plot["precison_with_unknown"]
		s = df_to_plot["sensitivity"]
		ax1.scatter(p, s, s=20, c="k", marker=marker_kmer[i], label=str(kmer)+"-mer",edgecolors='none',alpha=alpha )
	# plot legend after
	ax1.legend(loc=loc_label)

ax1.set_xlabel(x_label)
ax1.set_ylabel(y_label)

if manual_scale:
	ax1.set_xlim(xlim)
	ax1.set_ylim(ylim)


fig1.tight_layout()
fig1.subplots_adjust(top=0.92)

if save_output:
	fig1.savefig(filename_output, dpi=300)

else:
	fig1.show()


