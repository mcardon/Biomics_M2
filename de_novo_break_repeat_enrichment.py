# -*- coding: utf-8 -*-
# arg1 : file of filenames : csv files with positions of breaks (absolute path)
# arg2 : shustring_file.out or csv files with repeat data
# arg3 : type of repeat file : shustring or csv



################################ IMPORT ################################################################################################

#from sequana.lazy import pylab
from sequana.lazy import pandas as pd
#from sequana.lazy import numpy as np
#from random import shuffle
import sys
import os
from scipy.stats import chisquare

################################ PARAMETERS ##############################################################################################

fof_breaks_csv = str(sys.argv[1])
file_shustring = str(sys.argv[2])
type_repeat_file = str(sys.argv[3])

if type_repeat_file not in ["shustring","csv"]:
	print("arg 3 error : type of file not available : possible values : shustring csv")

threshold = 1000

################################ FUNCTIONS ##############################################################################################

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0



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


def rotate_breaks(df):
	# rotate the break if needed
	# rotate whole break
	ind_rotate = df[((df["begin"] % len_genome) != df["begin"]) & ((df["end"] % len_genome) != df["end"])].index.tolist()
	for i in ind_rotate:
		df.loc[i,("begin")] = df.loc[i,("begin")] % len_genome
		df.loc[i,("end")] = df.loc[i,("end")] % len_genome

	# if break bypass the start : cut in two parts
	ind_rotate = df[((df["begin"] % len_genome) == df["begin"]) & ((df["end"] % len_genome) != df["end"])].index.tolist()
	for i in ind_rotate:
		df.loc[i,("begin")] = df.loc[i,("begin")] % len_genome
		df.loc[i,("end")] = len_genome
		df.append(df.iloc[i,],ingnore_index=True)
		df.loc[df.shape[0],("begin")] = 0
		df.loc[df.shape[0],("end")] = df[i,("end")] % len_genome

	# sort by position and update index
	df.sort_values(by="begin",inplace=True)
	df.index = range(df.shape[0])


def find_overlap(start_break, end_break, start_rep, end_rep):
	start_overlap, end_overlap = max(start_break, start_rep), min(end_break, end_rep)
	if start_overlap > end_overlap:
		start_overlap, end_overlap = None, None
	return start_overlap, end_overlap


def find_nb_overlap_repeats(start_break, end_break, df_repeats):
	"""
	df_repeats must be sorted by begin position !!
	"""
	rep_i = 0
	nb_repeat_overlap = 0
	while rep_i < df_repeats.shape[0]:
		start_rep = df_repeats.loc[rep_i,"begin"]
		end_rep = df_repeats.loc[rep_i,"end"]

		# check if repeat can overlap the break
		if not (start_break > end_rep) :
			# check if start inside repeat
			# print("repeat : %d , %d" %(start_rep, end_rep))
			start_break_inside = (start_break >=  start_rep) & (start_break <= end_rep)
			end_break_inside = (end_break >=  start_rep) & (end_break <= end_rep)
			start_rep_inside = (start_rep >=  start_break) & (start_rep <= end_break)
			end_rep_inside = (end_rep >=  start_break) & (end_rep <= end_break)
			if start_break_inside | end_break_inside :
				# find overlap
				start_overlap, end_overlap = find_overlap(start_break, end_break, start_rep, end_rep)
				# print("overlap : %d , %d" %(start_overlap, end_overlap))
				nb_repeat_overlap += (end_overlap - start_overlap)
			elif start_rep_inside | end_rep_inside :
				# find overlap
				start_overlap, end_overlap = find_overlap(start_break, end_break, start_rep, end_rep)
				# print("overlap : %d , %d" %(start_overlap, end_overlap))
				nb_repeat_overlap += (end_overlap - start_overlap)
		# all repeats after this cannot overlap
		if end_break < start_rep:
			break
		rep_i +=1

	return nb_repeat_overlap

################################ INPUT DATA ##############################################################################################

print("Breakpoints files")

# list of all files
f = open(fof_breaks_csv, 'r')
list_files_breaks = []
list_df_breaks_files = []
for path in f.readlines():
	f_path = path.split('\n')[0]
	list_files_breaks.append(f_path)

	# save short name (for plots)
	short_name = f_path.split("/")[-1]
	list_df_breaks_files.append(short_name.replace(".csv",""))
f.close()


# load all breaks files and set header
list_df_breaks = []
list_df_breaks_names = []
for i in range(len(list_files_breaks)):
	f_path = list_files_breaks[i]
	if is_non_zero_file(f_path):
		# import df
		df_breaks = pd.read_csv(f_path,header=None,sep=",")
		df_breaks.columns = ["begin","end"]
		# add filename in df
		df_breaks["filename"] = [list_df_breaks_files[i]]*df_breaks.shape[0]
		# save in list
		list_df_breaks.append(df_breaks)
		list_df_breaks_names.append(list_df_breaks_files[i])
	else:
		print("%s : empty file, not imported" % list_df_breaks_files[i])



if type_repeat_file == "shustring":
	print("Shustring file")
	# load shustring file with pandas, and re-read to set header and find genome length
	df_shustring = pd.read_csv(file_shustring,comment='#',header=None,sep="\s+")

	f = open(file_shustring, 'r')
	l = f.readline()
	# first line contains genome length
	len_genome = int(l.split("\t")[1])
	# find header
	while(l):
		if (l[0] == '#') & (l[1] != '>'):
			header_df = l.replace('\n','').replace('#','').split('\t')
			break
		l = f.readline()
	f.close()
	df_shustring.columns = header_df

	print("Find repeats location and save as csv file")
	step_repeat_seq, be_repeats_concat = step_repeat_len(df_shustring, threshold)
	df_repeats = pd.DataFrame(be_repeats_concat)
	df_repeats.columns = ["begin","end"]
	path_outfile = file_shustring + "_df_begin_end.csv"
	# if csv file doesn't exist, save it
	if not os.path.isfile(path_outfile):
		f = open(path_outfile, 'a')
		f.write("# %d\n" %len_genome)
		df_repeats.to_csv(f,mode = 'a', index=None)
		f.close()

elif type_repeat_file == "csv":
	print("Repeat csv file")
	df_repeats = pd.read_csv(file_shustring,comment='#')

	f = open(file_shustring, 'r')
	l = f.readline()
	# first line contains genome length
	len_genome = int(l.replace("\n","").split(" ")[1])
	f.close()


################################ EXECUTE ##############################################################################################

# expected frequency of repeats in genome
freq_exp_inside_repeat = 0
for rep_i in range(df_repeats.shape[0]):
	freq_exp_inside_repeat += (df_repeats.loc[rep_i,"end"] - df_repeats.loc[rep_i,"begin"])
freq_exp_outside_repeat = len_genome - freq_exp_inside_repeat
# normalisation
freq_exp_inside_repeat = freq_exp_inside_repeat / float(len_genome)
freq_exp_outside_repeat = freq_exp_outside_repeat / float(len_genome)

# observed frequency of repeats between contigs (breaks)
for df_i in range(len(list_df_breaks)):
	df = list_df_breaks[df_i]
	df["nb_outside_repeat"] = [0]*df.shape[0]
	df["nb_inside_repeat"] = [0]*df.shape[0]

	# rotate breaks if needed
	rotate_breaks(df)
	
	for break_i in range(df.shape[0]):
		start_break = df.loc[break_i,"begin"]
		end_break = df.loc[break_i,"end"]

		nb_repeat_overlap = find_nb_overlap_repeats(start_break, end_break, df_repeats)

		df.loc[break_i,"nb_inside_repeat"] = nb_repeat_overlap
		df.loc[break_i,"nb_outside_repeat"] = (end_break - start_break) - nb_repeat_overlap

	# results repeats
	nb_obs_inside_repeat = df["nb_inside_repeat"].sum()
	nb_obs_outside_repeat = df["nb_outside_repeat"].sum()
	total = nb_obs_inside_repeat + nb_obs_outside_repeat
	# normalisation
	# nb_obs_inside_repeat = nb_obs_inside_repeat / float(total)
	# nb_obs_outside_repeat = nb_obs_outside_repeat / float(total)

	result_chi = chisquare(f_obs=[nb_obs_inside_repeat,nb_obs_outside_repeat],f_exp=[round(freq_exp_inside_repeat*total), round(freq_exp_outside_repeat*total)])
	print("\n" + list_df_breaks_names[df_i])
	print("Expected (repeat / non repeat) : ", [round(freq_exp_inside_repeat*total), round(freq_exp_outside_repeat*total)])
	print("Observed (repeat / non repeat) : ", [nb_obs_inside_repeat,nb_obs_outside_repeat])
	print(result_chi)



























