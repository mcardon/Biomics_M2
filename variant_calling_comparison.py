# -*- coding: utf-8 -*-
# arg1 : file of filenames : vcf files (from smrt : names will look better if : 'output_NN.vcf' where NN is the name of analysis)
# arg2 : file illumina (csv format)
# arg3 : shustring_file.out
# arg4 : filename for figure (if 'show', don't save but show instead)

# Only one file is allowed for illumina
# shustring file is the output of shustring tool 

################################ IMPORT ################################################################################################

from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
import collections
import sys

################################ PARAMETERS ##############################################################################################


file_of_filenames = str(sys.argv[1])
file_illumina = str(sys.argv[2])
file_shustring = str(sys.argv[3])
file_fig = str(sys.argv[4])

# threshold for repeat length
threshold = 100

# name of result file (withut extension)
filename_all_results = "2017_06_01_Results_variants_all"
filename_summary = "2017_06_01_Results_variants_summary"

# step to plot genome
step = 200000

# plot params
fontsize_x = 14
markersize = 6
hspace_subplots = 1.2
xtick_too_close = 15000

colormap = 'nipy_spectral_r'
# colormap = 'gist_rainbow_r'
# custom colors
custom_colors = True
# colors = ['k','r','b','g','b','c','k']
# colors =  ["#ff2d00","#ff2d00","#ff2d00","#ff2d00","#ccf900","#ccf900","#ccf900",\
# "#00ba00","#00ba00","#00ba00","#00a4bb","#00a4bb","#00a4bb",\
# "#0000b9","#0000b9","#0000b9","#000000","#000000","#000000"]

colors =  ["#ff0000","#ff0000","#005a9a","#00ba00","#00a4bb","#0000b9","#000000"]


################################ FUNCTIONS ################################################################################################


def step_repeat_len(df_shustring, threshold):
	"""
	Convert repetitivity curve into step curve with repeat len
	returns begining and end as list of tuples
	"""

	nb_row = df_shustring.shape[0]
	i = 0
	step_repeat_seq = []
	be_repeats = []
	e = 0

	# use list because much faster
	list_len_shus = list(df_shustring.iloc[:,1])

	while(i < nb_row):
		# begining of repeat
		if (list_len_shus[i] > threshold):
			b = i
			#print(b)
			# add 0 for non repeat area
			step_repeat_seq = step_repeat_seq + [0 for j in range(e,b,1)]
			# compute new end of repeat
			len_repeat = list_len_shus[i]
			e = b + len_repeat
			# add len of repeat for all repeat area
			step_repeat_seq = step_repeat_seq + [len_repeat for j in range(b,e,1)]
			# save (b,e)
			be_repeats.append((b,e))
			# update i
			i = e-1
		i +=1

	step_repeat_seq = step_repeat_seq + [0 for j in range(e,i,1)]

	# merge repeats that are fragmented
	prev_tup = be_repeats[0]
	b = prev_tup[0]
	be_repeats_concat = []
	for tup in be_repeats[1:len(be_repeats)]:
		if tup[0] == prev_tup[1]:
			# concat
			e = tup[1]
			prev_tup = tup
		else:
			# real end of repeat : append result and update b, e
			e = prev_tup[1]
			be_repeats_concat.append((b,e))
			prev_tup = tup
			b = prev_tup[0]

	return step_repeat_seq, be_repeats_concat


def convert_illumina_indel(small_seq, large_seq, indel):
	"""
	Converts insertion notation from variant calling sequana to match notation smrt pacbio
	"""
	small_seq = list(small_seq[1:(len(small_seq)-1)])
	large_seq = list(large_seq[1:(len(large_seq)-1)])
	while(small_seq):
		if small_seq[-1] == large_seq[-1]:
			small_seq.pop()
			large_seq.pop()

	if indel == 'insert':
		new_small_seq = '.'
		new_large_seq = 'I' + ''.join(large_seq)
	elif indel == 'delete':
		new_small_seq = 'D' + str(len(large_seq))
		new_large_seq = ''.join(large_seq)

	return new_small_seq, new_large_seq


def type_variants_illumina(df_illumina, column_CDS):
	# is it a mutation or indel ?
	type_var = list(df_illumina[column_CDS].values)
	mut_indel = []
	for cds in type_var:
		if ">" in cds:
			mut_indel.append('mut')
		else:
			if ("ins" in cds) & ("del" in cds):
				mut_indel.append('indel')
			elif "ins" in cds:
				mut_indel.append('ins')
			elif "del" in cds:
				mut_indel.append('del')
			else:
				mut_indel.append('unknown')
	return mut_indel

def type_variants_smrt(df_smrt, column_ALT):
	type_var = list(df_smrt[column_ALT].values)
	mut_indel = []
	for cds in type_var:
		if ("I" in cds) & ("D" in cds):
			mut_indel.append('indel')
		elif "I" in cds:
			mut_indel.append('ins')
		elif"D" in cds:
			mut_indel.append('del')
		else:
			mut_indel.append('mut')
	return mut_indel

def compare_results(df_result ,analysis_names,what_to_compare):
	pacbio_identical = [True]*df_result.shape[0]
	to_compare = []
	if what_to_compare == 'pacbio':
		start = 1
	else:
		start = 0
	for analysis in analysis_names[start:len(analysis_names)]:
		if list(to_compare):
			#print(analysis)
			identical = to_compare == df_result[analysis]
			pacbio_identical = pacbio_identical & identical
			to_compare = df_result[analysis]
		else:
			to_compare = df_result[analysis]
	return pacbio_identical

def summary_results(df_analysis, name):
	print(name, df_analysis.shape)
	res_to_append = []
	# analysis name 0
	res_to_append.append(name)
	#number of variants found 1
	res_to_append.append(df_analysis.shape[0])
	# number outside repeat 2
	res_to_append.append(df_analysis[df_analysis['size_rep'] == 0].shape[0])
	# number inside repeat 3
	df_in_rep = df_analysis[df_analysis['size_rep'] != 0]
	res_to_append.append(df_in_rep.shape[0])
	# mean and median rep size 4 5
	res_to_append.append(df_in_rep['size_rep'].mean())
	res_to_append.append(df_in_rep['size_rep'].median())
	# min and max rep size 6 7
	res_to_append.append(df_in_rep['size_rep'].min())
	res_to_append.append(df_in_rep['size_rep'].max())
	# type of variant mut ins del 8 9 10
	res_to_append.append(df_analysis[df_analysis['mut_indel'] == 'mut'].shape[0] )
	res_to_append.append(df_analysis[df_analysis['mut_indel'] == 'ins'].shape[0] )
	res_to_append.append(df_analysis[df_analysis['mut_indel'] == 'del'].shape[0] )
	not_sure = df_analysis.shape[0] - sum(res_to_append[-3:])
	res_to_append.append(not_sure)

	return res_to_append


################################ FUNCTIONS FOR PLOTS #########################################################################################

def find_begin_end_interval(be_repeats_concat, i, i_end):
	"""
	Returns list of repeats in given interval
	"""

	# find begin end of repeats in interval
	be_repeats_concat_i = [j for j in be_repeats_concat if (((j[0] > i) & (j[0] < i_end)) | ((j[1] > i) & (j[1] < i_end)))]

	# check first and last element, trim edges if necessary
	if len(be_repeats_concat_i) > 1 :
		r_first = be_repeats_concat_i.pop(0)
		r_last = be_repeats_concat_i.pop()
		be_repeats_concat_i = [( max(r_first[0], i) ,r_first[1])] + be_repeats_concat_i + [(r_last[0], min(r_last[1], i_end ))]
	elif len(be_repeats_concat_i) == 1:
		r_last = be_repeats_concat_i.pop()
		be_repeats_concat_i = [( max(r_last[0],i) , min(r_last[1],i_end) )]

	return be_repeats_concat_i


def subplot_variant_position(df_result, i, gen_pos, axarr, analysis_names, y_pos, y_col, be_repeats_concat):
	ax = axarr[i]
	pos = gen_pos[i]
	ax.set_xlim(pos)
	ax.set_ylim([y_pos[0],y_pos[-1]])
	df_to_plot = df_result[(df_result['position'] >= pos[0]) & (df_result['position'] <= pos[1]) ]
	my_xticks = [pos[0]] + list(df_to_plot['position']) + [pos[1]]

	# merge ticks that are too close
	prev_t = pos[0]
	to_append = str(prev_t)
	new_my_xticks = []
	to_pop = []
	for j in range(1,len(my_xticks[1:len(my_xticks)]),1):
		t = my_xticks[j]
		if abs(t - prev_t) < xtick_too_close:
			to_append = to_append + "\n " + str(t)
			prev_t = t
			to_pop.append(j)
		else:
			new_my_xticks.append(to_append)
			prev_t = t
			to_append = str(prev_t)
	new_my_xticks.append(to_append)

	ticks_pos = [my_xticks[j] for j in range(len(my_xticks)) if j not in to_pop]
			

	ax.set_xticks(ticks_pos)
	ax.set_xticklabels(new_my_xticks, fontsize=fontsize_x)
	for j in range(len(analysis_names)):
		analysis = analysis_names[j]
		df = df_to_plot.dropna(axis=0, how='all', subset=[analysis])
		if analysis == 'Illumina':
			ax.plot(df['position'], [y_pos[j+1]]*df.shape[0], 'ks',markersize=markersize, label=analysis)
		else:
			ax.plot(df['position'], [y_pos[j+1]]*df.shape[0], color=y_col[j], marker='o', linestyle='None', label=analysis )
	# highlight repeated regions

	be_repeats_concat_i = find_begin_end_interval(be_repeats_concat, gen_pos[0][0], gen_pos[0][1])
	for rep in be_repeats_concat:
		ax.axvspan(rep[0], rep[1], alpha=0.5, color='orange')



################################ INPUT DATA ##############################################################################################

# list of all files
f = open(file_of_filenames, 'r')
list_files_smrt = []
list_df_smrt_names = []
for name in f.readlines():
	f_name = name.split('\n')[0]
	list_files_smrt.append(f_name)

	# save short name (for plots)
	short_name = f_name.split("/")[-1]
	list_df_smrt_names.append(short_name)
f.close()


# load all smrt files and set header
print("\nSMRT files")
list_df_smrt = []

for f_smrt in list_files_smrt:
	# import df
	df_smrt = pd.read_csv(f_smrt,comment='#',header=None,sep="\t")

	# find header
	f = open(f_smrt, 'r')
	l = f.readline()
	while(l):
		if (l[0] == '#') & (l[1] != '#'):
			header_df = l.replace('\n','').replace('#','').split('\t')
			break
		l = f.readline()
	f.close()
	df_smrt.columns = header_df

	# save in list
	list_df_smrt.append(df_smrt)


# load illumina file
print("\nIllumina file")
df_illumina = pd.read_csv(file_illumina,sep=",")

print("Shustring file")
# load shustring file and set header
df_shustring = pd.read_csv(file_shustring,comment='#',header=None,sep="\s+")
# find header
f = open(file_shustring, 'r')
l = f.readline()
# first line contains genome length
len_genome = int(l.split("\t")[1])

while(l):
	if (l[0] == '#') & (l[1] != '>'):
		header_df = l.replace('\n','').replace('#','').split('\t')
		break
	l = f.readline()
f.close()
df_shustring.columns = header_df


################################ EXECUTE : ALL RESULTS ##############################################################################################

print("find repeats location")
# find repeats location
step_repeat_seq, be_repeats_concat = step_repeat_len(df_shustring, threshold)

# merge data from Illumina and Pacbio smrt
## annotate before merge
print("Merge data")
analysis_names = []

# Illumina
df_illumina['data_type'] = ['Illumina']*df_illumina.shape[0]
ref_list = list(df_illumina['reference'].values)
alt_list = list(df_illumina['alternative'].values)
tri_del = [len(ref_list[i]) > len(alt_list[i]) for i in range(len(alt_list))]
# add 1 to position in variant calling sequana, to match position of smrtlink
df_illumina.loc[tri_del,'position'] =  df_illumina.loc[tri_del,'position']  +1

# add column with type of mutation
mut_indel = type_variants_illumina(df_illumina, 'CDS_position')
df_illumina['mut_indel'] = mut_indel

# append to list of results
df_to_merge = df_illumina[['position','reference','alternative','freebayes_score','mut_indel','data_type']]
df_to_merge.columns = ['position','reference','alternative','score','mut_indel','data_type']
list_df_to_merge = [df_to_merge]
analysis_names.append('Illumina')

#Pacbio
for i in range(len(list_df_smrt_names)):
	df_smrt = list_df_smrt[i]
	name = list_df_smrt_names[i]
	name = name.replace('output_','').replace('.vcf','')
	analysis_names.append('Pacbio_' + name)
	# annotate
	df_smrt['data_type'] = ['Pacbio_' + name]*df_smrt.shape[0]

	# add column with type of mutation
	mut_indel = type_variants_smrt(df_smrt, 'ALT')
	df_smrt['mut_indel'] = mut_indel

	df_to_merge = df_smrt[['POS','REF','ALT','QUAL','mut_indel','data_type']]
	df_to_merge.columns = ['position','reference','alternative','score','mut_indel','data_type']

	list_df_to_merge.append(df_to_merge)


concat_variant = pd.concat(list_df_to_merge)
concat_variant.sort_values(by='position',axis=0,inplace=True)


print("Create results table")
# format columns names
score_names = [ analysis + '_score' for analysis in analysis_names]
col_names = []
for i in range(len(analysis_names)):
	col_names.extend([analysis_names[i], score_names[i]])

res_variant_names = ['position'] + col_names + ['mut_indel']
res_variant = []

for i in set(concat_variant['position']):
	df = concat_variant[concat_variant['position'] == i]
	# position
	res = [i]
	# analysis and score
	var_found = set(df['data_type'])
	for j in analysis_names:
		if j in var_found:
			ref = df.loc[df['data_type'] == j, 'reference'].values[0]
			alt = df.loc[df['data_type'] == j, 'alternative'].values[0]
			score = df.loc[df['data_type'] == j, 'score'].values[0]
			# indel in illumina are different format from pacbio, change to compare
			if (j == "Illumina") & (len(ref) != len(alt)):
				if len(ref) > len(alt):
					alt, ref = convert_illumina_indel(alt, ref, 'delete')
				else:
					ref, alt = convert_illumina_indel(ref, alt, 'insert')

			res.append(ref + '>' + alt)
			res.append(score)
		else:
			res.append(None)
			res.append(None)
	# type of mutation : should be unique, but consider it's not, just in case
	var_type = set(df['mut_indel'])
	res.append('_'.join(var_type))

	# append result
	res_variant.append(res)


df_result = pd.DataFrame(res_variant)
df_result.columns = res_variant_names # XXX change names
df_result.sort_values(by='position', axis=0, inplace=True)


df_result['identical_pacbio'] = compare_results(df_result ,analysis_names, 'pacbio')
df_result['identical_all'] = compare_results(df_result ,analysis_names, 'all')

# for each position, check if inside repeat
be_deque = collections.deque(be_repeats_concat)
prev_rep_end = 0
next_rep = be_deque.popleft()
is_inside_repeat = []
size_rep = []
for pos in list(df_result['position']):
	while not (pos >= prev_rep_end) & (pos <= next_rep[1]):
		# update deque rep and previous position
		prev_rep_end = next_rep[1]
		if be_deque:
			next_rep = be_deque.popleft()
		else:
			#deque is empty, set ed of genome
			next_rep = (len_genome,len_genome)

	# check if inside rep
	if pos >= next_rep[0]:
		is_inside_repeat.append(True)
		size_rep.append(next_rep[1] - next_rep[0])
	else:
		is_inside_repeat.append(False)
		size_rep.append(0)
# append result to df
df_result['is_inside_repeat_' +str(threshold)] = is_inside_repeat
df_result['size_rep'] = size_rep

print("Save table of all results")
df_result.to_csv(filename_all_results + ".csv")


################################ EXECUTE : SUMMARY ##############################################################################################

print("Create summary")
header_df = ['analysis','total_variants','variants_outside_repeat', 'variant_inside_repeat', 'mean_repeat_size', 'median_repeat_size', 'min_repeat_size', 'max_repeat_size', 'mut','ins','del','not_sure']
result = []
for analysis in analysis_names:
	df_analysis = df_result[df_result[analysis].notnull()]
	res_to_append = summary_results(df_analysis, analysis)

	# append all
	result.append(res_to_append)

# All pacbio analysis (union)
Pacbio_cols = [col for col in df_result.columns if 'Pacbio' in col]
df_analysis = df_result.dropna(axis=0,how='all',subset=Pacbio_cols)
res_to_append = summary_results(df_analysis, 'Union_pacbio')
result.append(res_to_append)

# Only variants found in all pacbio (intersection)
df_analysis = df_result[df_result['identical_pacbio'] == True]
res_to_append = summary_results(df_analysis, 'Intersection_pacbio')
result.append(res_to_append)

# Only variants found in all analysis (intersection all analysis)
df_analysis = df_result[df_result['identical_all'] == True]
res_to_append = summary_results(df_analysis, 'Intersection_all')
result.append(res_to_append)

# Found in pacbio but not illumina
df_analysis = df_result[df_result['Illumina'].isnull()]
res_to_append = summary_results(df_analysis, 'Only_pacbio')
result.append(res_to_append)

# Found in illumina but not pacbio
dont_keep = df_result.dropna(axis=0,how='all',subset=Pacbio_cols)
df_analysis = df_result[~df_result.index.isin(dont_keep.index)]
res_to_append = summary_results(df_analysis, 'Only_Illumina')
result.append(res_to_append)

summary_variants = pd.DataFrame(result)
summary_variants.columns = header_df

print("Save summary")
summary_variants.to_csv(filename_summary + ".csv")



################################ PLOTS ##############################################################################################


print("Create plots")
#colors = ['m','r','y','g','b','c','k']

cmap = pylab.cm.get_cmap(colormap)
#colors = [cmap(i) for i in np.linspace(0,1,len(list_analysis))]

# positions of genome
gen_pos = [[i, i+step-1] for i in range(0,len_genome,step)]
y_pos = list(np.linspace(0,1,len(analysis_names)+2))

if custom_colors:
	y_col = [colors[i] for i in range(len(analysis_names))]
else:
	y_col = [cmap(i) for i in np.linspace(0,1,len(analysis_names)) ]

pylab.close('all')

# create figure
fig, axarr = pylab.subplots(len(gen_pos),1, figsize=(int(step/20000), int(len(gen_pos))*1.1))
for i in range(len(gen_pos)):
	subplot_variant_position(df_result, i, gen_pos, axarr, analysis_names, y_pos, y_col, be_repeats_concat)

# add grey at the end (no genome)
ax = axarr[-1]
ax.axvspan(len_genome,gen_pos[-1][1], alpha=0.5, color='k')

#fig.subplots_adjust(bottom=0.2)
#fig.tight_layout()
pylab.subplots_adjust(hspace = hspace_subplots)
pylab.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#pylab.legend(loc="lower left")
if file_fig != "show":
	pylab.savefig(file_fig)
else:
	pylab.show()
pylab.close('all')
#pylab.show()








