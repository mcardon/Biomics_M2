# -*- coding: utf-8 -*-
# arg1 : what to do : possible choices : covplot hist2d boxplot all
# arg2 : filename for plot (if "show", dont save, show instead)

################################ IMPORT ################################################################################################

import pandas as pd
import pylab
import sys
from biokit.viz import hist2d

################################ PARAMETERS ##############################################################################################

# sys parameters
if len(sys.argv) > 0:
	# what to do
	what_to_do = str(sys.argv[1])

	do_coverage_plot = False
	do_coverage_hist2d = False
	do_boxplots = False

	if what_to_do == 'covplot':
		do_coverage_plot = True
	elif what_to_do == 'hist2d':
		do_coverage_hist2d = True
	elif what_to_do == 'boxplot':
		do_boxplots = True
	elif what_to_do == 'all':
		do_coverage_plot = True
		do_coverage_hist2d = True
		do_boxplots = True
	else:
		print('Incorrect arg : %s \nPossible choices : covplot hist2d boxplot all' %what_to_do)

	# save or show
	save_plots = str(sys.argv[2])
	if save_plots == 'save':
		save_plots = True
	elif save_plots == 'show':
		save_plots = False
	else:
		print('Incorrect arg : %s \nPossible choices : save show' %save_plots)


# step to plot genome
step = 200000

# lenght threshold for repeated sequence
threshold = 100

# prefix name plots
dir_save = "/home/mcardon/Mel/Rapport/Figures/"
figname_cov     = dir_save + "Coverage_plot_rep" +str(threshold)
figname_hist2D  = dir_save + "Hist2D_cov_repeat_rep" +str(threshold)
figname_boxplot = dir_save + "Boxplot_cov_repeat_len_rep" +str(threshold)

# display parameters
fontsize = 12


################################ LOAD DATA ################################################################################################

# ref genome shustring (output without header)
genom_ref = pd.read_csv("/home/mcardon/Mel/Datas/Ref_genome_NC_002929/NC_002929_shustring_nohead.out",header=None,sep="\s+")

# bed files wit coverage data
cov_pacbio_blasr = pd.read_csv("/home/mcardon/Mel/Datas/Ref_genome_NC_002929/cov_pacbio_blasr_mapq100.bed",header=None,sep="\s+")
cov_pacbio_bwa = pd.read_csv("/home/mcardon/Mel/Datas/Ref_genome_NC_002929/cov_pacbio_bwa_mem_mapq30.bed",header=None,sep="\s+")
cov_illumina = pd.read_csv("/home/mcardon/Mel/Datas/Ref_genome_NC_002929/cov_illumina_mapq30.bed",header=None,sep="\s+")

print("Data loaded")

################################ FUNCTIONS ################################################################################################


def find_repeat_location(genom_ref, threshold):
	"""
	Returns list of tuple [(begin_repeat, end_repeat)]

	note = not optimised. Would be better to use a stack
	if empty stack, add location
	if already a location in stack, pop it (this is the begining) and save (b,e) in list
	"""
	repeat_seq = [1 if (i > threshold) else 0 for i in genom_ref[1] ]

	out_rep = True
	repeat_pos = []
	for i in range(len(repeat_seq)):
		# begining of repeat
		if out_rep & repeat_seq[i]:
			b = i
			out_rep  = False

		# end of repeat
		if (not out_rep) & (not repeat_seq[i]):
			out_rep = True
			e = i
			repeat_pos.append((b,e))

	return repeat_pos


def find_begin_end_interval(repeat_pos, i, i_end):
	"""
	Returns list of repeats in given interval
	"""

	# find begin end of repeats in interval
	repeat_pos_i = [j for j in repeat_pos if (((j[0] > i) & (j[0] < i_end)) | ((j[1] > i) & (j[1] < i_end)))]

	# check first and last element, trim edges if necessary
	if len(repeat_pos_i) > 1 :
		r_first = repeat_pos_i.pop(0)
		r_last = repeat_pos_i.pop()
		repeat_pos_i = [( max(r_first[0], i) ,r_first[1])] + repeat_pos_i + [(r_last[0], min(r_last[1], i_end ))]
	elif len(repeat_pos_i) == 1:
		r_last = repeat_pos_i.pop()
		repeat_pos_i = [( max(r_last[0],i) , min(r_last[1],i_end) )]

	return repeat_pos_i


def df_filter_repetitivity(cov_df, genom_ref,threshold):
	"""
	Returns 2 df with lw repetitivity and high repetitivity data
	"""
	# complete df
	nb_row = genom_ref.shape[0]
	df = pd.DataFrame(index=genom_ref[0], columns= ['cov','repetitivity'])
	df['cov'] = cov_df.iloc[0:nb_row,2]
	df['repetitivity'] = genom_ref[1]
	
	df_low_rep = df[df.repetitivity < threshold]
	df_high_rep = df[df.repetitivity >= threshold]
	return df_low_rep, df_high_rep



def plot_hist2D_cov(df_low_rep, df_high_rep, name_fig, title_fig, save_plots):
	"""
	"""
	### histogram 2D
	# low repetitivity
	h = hist2d.Hist2D(df_low_rep)
	h.plot(bins=[threshold,10], cmap='gnuplot2_r',fontsize=fontsize, xlabel='Coverage', ylabel='Repetitivity')
	pylab.title('%s\n Low repetitivity genome' %title_fig)
	if save_plots:
		pylab.savefig("%s_%s_low_repetitivity.png" %(figname_hist2D,name_fig))
	else:
		pylab.show()

	# high repetitivity
	h = hist2d.Hist2D(df_high_rep)
	h.plot(bins=[30,300], cmap='gnuplot2_r', fontsize=fontsize, xlabel='Coverage', ylabel='Log Repetitivity')
	pylab.yscale('log')
	pylab.title('%s\n High repetitivity genome' %title_fig)
	if save_plots:
		pylab.savefig("%s_%s_high_repetitivity.png" %(figname_hist2D,name_fig))
	else:
		pylab.show()


def step_repeat_len(genom_ref, threshold):
	"""
	Convert repetitivity curve into step curve with repeat len
	"""

	nb_row = genom_ref.shape[0]
	i = 0
	step_repeat_seq = []
	e = 0

	# use list because much faster
	list_len_shus = list(genom_ref[1])#[i]

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
			# update i
			i = e-1
		i +=1

	step_repeat_seq = step_repeat_seq + [0 for j in range(e,i,1)]

	return step_repeat_seq

def boxplot_cov_repeat_len(genom_ref, levels_repeat, cov_data, name_fig, title, save_plots):
	"""
	"""
	nb_row = genom_ref.shape[0]
	max_repeat = max(step_repeat_seq)

	cov = cov_data.iloc[0:nb_row,:] 
	data_boxplot = [cov[ (genom_ref.step_repeat >= i[0]) & (genom_ref.step_repeat <= i[1]) ][2] for i in levels_repeat]
	x_labels = [str(levels_repeat[i]) + "\nN = " + str(len(data_boxplot[i])) for i in range(len(levels_repeat))]

	# Boxplot coverage
	pylab.close('all')
	pylab.figure(figsize=(len(levels_repeat)*1.2, 5))
	pylab.boxplot(data_boxplot)
	pylab.xticks(range(1,len(levels_repeat)+1,1),x_labels)
	pylab.xlabel("Repeat length")
	pylab.ylabel("Coverage")
	pylab.title("Coverage vs Repeat length\n%s" % title)
	pylab.tight_layout()
	if save_plots:
		pylab.savefig("%s_%s.png" % (figname_boxplot, name_fig))
	else:
		pylab.show()
	pylab.close()


################################ EXECUTE ################################################################################################

# genome size
gen_size = genom_ref.shape[0]

repeat_pos = find_repeat_location(genom_ref, threshold)

step_repeat_seq = step_repeat_len(genom_ref, threshold)
genom_ref['step_repeat'] = step_repeat_seq

################################ PLOTS ################################################################################################

########################### Coverage plots with hightlight on repeated sequences
if do_coverage_plot:
	for i in range(0,gen_size,step):
		i_end = min(i+step, gen_size -1)

		# find begin end of repeats in interval
		repeat_pos_i = find_begin_end_interval(repeat_pos, i, i_end)

		pylab.close('all')
		# create figure
		fig, axarr = pylab.subplots(4,1, sharex=True, figsize=(step*0.00008, 5))

		#ax = pylab.subplot(411, sharex=True)
		axarr[0].set_title("Coverage Illumina with bwa mem (mapq >30)")
		axarr[0].plot(cov_illumina[2][i:i_end],lw=0.7)
		
		#axarr = pylab.subplot(412, sharex=True)
		axarr[1].plot(cov_pacbio_bwa[2][i:i_end],lw=0.3)
		axarr[1].set_title("Coverage Pacbio with bwa mem (mapq > 30)")

		#axarr = pylab.subplot(413, sharex=True)
		axarr[2].plot(cov_pacbio_blasr[2][i:i_end],lw=0.3)
		axarr[2].set_title("Coverage Pacbio with blasr (mapq > 100)")

		#axarr = pylab.subplot(414, sharex=True)
		axarr[3].set_title("Repeat index")
		axarr[3].plot(genom_ref[1][i:i_end])
		#axarr[3].yscale('log')

		for rep in repeat_pos_i:
			for j in range(len(axarr)):
				axarr[j].axvspan(rep[0], rep[1], alpha=0.5, color='orange')

		pylab.tight_layout()
		if save_plots:
			pylab.savefig("%s_%d_%d.png" %(figname_cov,i,i_end))
		else:
			pylab.show()
		pylab.clf()
		pylab.close('all')



########################### hist2d : correlation between shustring length and coverage
if do_coverage_hist2d:
	pylab.close('all')

	##### Illumina data
	df_low_rep, df_high_rep = df_filter_repetitivity(cov_illumina, genom_ref,threshold)
	plot_hist2D_cov(df_low_rep, df_high_rep, "Illumina", "Illumina with bwa mem (mapq > 30)", save_plots)

	##### Pacbio bwa data
	df_low_rep, df_high_rep = df_filter_repetitivity(cov_pacbio_bwa, genom_ref,threshold)
	plot_hist2D_cov(df_low_rep, df_high_rep, "Pacbio_bwa", "Pacbio with bwa mem (mapq > 30)", save_plots)

	##### Pacbio blasr data
	df_low_rep, df_high_rep = df_filter_repetitivity(cov_pacbio_blasr, genom_ref,threshold)
	plot_hist2D_cov(df_low_rep, df_high_rep, "Pacbio_blasr", "Pacbio with blasr (mapq > 100)", save_plots)


########################### boxplot : coverage for repeat seq of different length

if do_boxplots:
	max_repeat = max(step_repeat_seq)	
	levels_repeat = [(0,49)] + [(i, i+99) for i in range(50, 250,100)] + [(i, i+249) for i in range(250, 1000,250)] + [(i, i+999) for i in range(1000,max_repeat,1000)]

	# Illumina data
	title = "Illumina with bwa mem (mapq > 30)"
	name_fig = "Illumina"
	boxplot_cov_repeat_len(genom_ref, levels_repeat, cov_illumina, name_fig, title, save_plots)

	# Pacbio data with bwa
	title = "Pacbio with bwa mem (mapq > 30)"
	name_fig = "Pacbio_bwa"
	boxplot_cov_repeat_len(genom_ref, levels_repeat, cov_pacbio_bwa, name_fig, title, save_plots)

	# Pacbio data with blasr
	title = "Pacbio with blasr (mapq > 100)"
	name_fig = "Pacbio_blasr"
	boxplot_cov_repeat_len(genom_ref, levels_repeat, cov_pacbio_blasr, name_fig, title, save_plots)








