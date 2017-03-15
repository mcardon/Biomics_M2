# -*- coding: utf-8 -*-
# arg1 : results_variant_all.csv : output of variant_calling_comparison.py
# arg2 : true_variants.csv : csv file with true variant info : must have 2 columns : 'position' 'variant'
# arg3 : length of genome
# arg4 : filename for figure (if 'show', don't save but show instead)

# only columns with "Pacbio" or "Illumina" in column name are taken into account for analysis

################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sklearn.metrics import precision_recall_curve
import sys

################################ PARAMETERS ##############################################################################################


file_result = str(sys.argv[1])
file_variants = str(sys.argv[2])
len_genome = int(sys.argv[3])
file_fig = str(sys.argv[4])


################################ FUNCTIONS ################################################################################################


def compute_table_performance(analysis, df_results):
	"""
	return [TP, FP, FN, TN]
	TP and FP are lists of scores
	"""
	TP = []
	FP = []
	FN = 0
	TN = 0
	df_analysis = df_results.loc[:,["position", analysis, analysis + "_score"]].dropna()

	# check all positions found
	pos_found_set = set(df_analysis["position"])

	# True positives : check that the variant is the same as expected
	# if not it's a False positive
	TP_candidates = pos_found_set & real_pos_set
	for pos in TP_candidates:
		var_found = df_analysis[df_analysis["position"] == pos].loc[:,analysis].values[0]
		var_real = df_variants[df_variants["position"] == pos].loc[:,"real_variant"].values[0]
		score = df_analysis[df_analysis["position"] == pos].loc[:,analysis + "_score"].values[0]
		if var_found == var_real:
			TP.append(score)
		else:
			FP.append(score)
			print("%s : there is %d right position found with wrong mutation" % (analysis,len(FP)))
	# false positives : fond as variants, but not true
	FP_set = pos_found_set - real_pos_set
	FP.extend([df_analysis[df_analysis["position"] == pos].loc[:,analysis + "_score"].values[0] for pos in FP_set])
	# false negatives : variants not found in analysis
	FN += len(real_pos_set - pos_found_set)
	# true positive : all positions found nowhere
	TN += len_genome - sum([len(TP), len(FP), FN])

	return [TP, FP, FN, TN]




################################ INPUT DATA ################################################################################################

df_results = pd.read_csv(file_result,sep=",")
df_variants = pd.read_csv(file_variants,sep=",")

################################ FUNCTIONS ################################################################################################

"""
##### XXX create list of true variant : should be done outside for final version (this will be an input of script)
cols = list(df_variants.columns)
col_to_keep = [ col for col in cols if ( ('Illumina' in col) | ('Pacbio' in col) | ('position' in col) )]
df_variants = df_variants[col_to_keep]
list_analysis = [col for col in col_to_keep if not ( ('score' in col) | ('position' in col) )]
list_pacbio_analysis = [col for col in list_analysis if ('Pacbio' in col)]

# if score < threshold, remove variant from analysis
for analysis in list_analysis:
	below_thr = df_variants[analysis + "_score"] < threshold
	rows_below_thr = [i for i in range(len(below_thr)) if below_thr[i]]
	df_variants.loc[rows_below_thr,analysis] = None


# variant is true if found in 4 out of 5 pacbio analysis
index_keep = []
for analysis in list_pacbio_analysis:
	list_col_to_test = [col for col in list_pacbio_analysis if col != analysis]
	print(list_col_to_test)
	for col_test in list_col_to_test:
		df = df_variants.loc[:,list_col_to_test]
		# get index where 4 analyses give same result
		index_keep.extend(list(df[df.apply(pd.Series.nunique, axis=1) == 1].dropna().index))

# unique values in ascending order
index_keep = sorted(list(set(index_keep)))
# remove rows with false positive
df_variants = df_variants.loc[index_keep,:]
# add column with real variant
df_variants["real_variant"] = df_variants.loc[:,list_pacbio_analysis].mode(axis=1)
df_variants = df_variants.loc[:,["position","real_variant"]]
df_variants.to_csv("Real_variants.csv")
"""

real_pos_set = set(df_variants["position"])

# get analysis names
cols = list(df_results.columns)
list_analysis = [ col for col in cols if ( (('Illumina' in col) | ('Pacbio' in col)) & ('score' not in col) )]

# change scores to have between 0 and 1
# for illumina : normalise with max
df_results["Illumina_score"] = df_results["Illumina_score"] / float(max(df_results["Illumina_score"].dropna()))
list_pacbio_analysis = [col for col in list_analysis if ('Pacbio' in col)]
for analysis in list_pacbio_analysis:
	df_results[analysis + "_score"] = round(df_results[analysis + "_score"] /100. , 2)


cmap = pylab.cm.get_cmap('plasma_r')
colors = [cmap(i) for i in np.linspace(0,1,len(list_analysis))]

# get results for curves
for i in range(len(list_analysis)):
	analysis = list_analysis[i]
	res = compute_table_performance(analysis, df_results)
	print("\n%s" %analysis)
	# [TP, FP, FN, TN]
	# print(len(res[0]), len(res[1]), res[2], res[3] , sum([len(res[0]), len(res[1]), res[2], res[3]]))
	TP = res[0]
	FP = res[1]
	FN = [0]*res[2]
	TN = [0]*res[3]
	y_true = np.array([1]*len(TP) + [1]*len(FN) + [0]*len(FP) + [0]*len(TN))
	y_scores = np.array(TP + FN + FP + TN)
	precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
	# print( precision  )
	# print( recall)
	# print( thresholds)
	pylab.plot(recall, precision, color=colors[i],label=analysis)
	pylab.xlabel('Recall')
	pylab.ylabel('Precision')
	pylab.ylim([0.0, 1.05])
	pylab.xlim([0.0, 1.05])
	pylab.title('Precision-Recall')
	pylab.legend(loc="lower left")

if file_fig != "show":
	pylab.savefig(file_fig)
else:
	pylab.show()


