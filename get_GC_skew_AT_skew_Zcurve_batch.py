# -*- coding: utf-8 -*-
# !! only for single chromosome (bact genome) 
# arg1 : file of filenames
# arg2 : file of filenames already done
# arg3 : type of window : fixed adapt
# arg4 : window size : int (if fixed) or float (proportion of genome)
# arg5 : prefix for figures names (optionnal, deflaut "./") if show : do not save, show instead

#### Ajouter GC content in title

################################ IMPORT ################################################################################################
import pylab
from Bio import SeqIO
import numpy as np
import collections
import sys

################################ PARAMETERS ##############################################################################################


file_of_filenames = str(sys.argv[1])
file_of_filenames_done = str(sys.argv[2])
type_window = str(sys.argv[3])

if type_window == 'fixed':
	window = int(sys.argv[4])
elif type_window == 'adapt':
	window = 0
	frac_window = float(sys.argv[4])
else:
	print("Error incorrect arg 3 : type of window. Possible values : fixed adapt")

if len(sys.argv) > 5:
	prefix_figure = str(sys.argv[5])
else:
	prefix_figure = "./"


# filename_ref = "/home/mcardon/Mel/Datas/Ref_genome_NC/NC_002929_Bpertussis.fasta"
# window = 50000 

################################ FUNCTIONS ##############################################################################################

def read_genome(filename_ref):
	"""
	"""
	# read genome reference
	ref_fasta = SeqIO.parse(filename_ref, "fasta")
	sequence_reference = str(next(ref_fasta).seq).upper()
	len_genome = len(sequence_reference)
	return sequence_reference, len_genome


def init_sliding_window(sequence_reference, window):
	"""
	init sliding window

	return :
	slide_window : deque of n first nucleotides (size of window)
	seq_right : deque of the rest of the genome, with the n-1 nucleotides at the end (circular DNA)
	"""

	# split sequence_reference in two : sliding window and seq_right
	# for circular genome : add begining at the end of genome
	slide_window = collections.deque(sequence_reference[0:window])
	seq_right = collections.deque( sequence_reference[window:len(sequence_reference)] + sequence_reference[0:(window-1)] )

	if len(sequence_reference) < window:
		print("Error : ref genome shorter than window")
		return 0
	else:
		return slide_window, seq_right


def init_list_results(len_genome):
	"""
	init empty np.array
	"""
	# init IJ content and IJ skew
	IJ_content_res = np.empty((1,len_genome))
	IJ_content_res[:] = np.NAN
	IJ_skew_res = np.empty((1,len_genome))
	IJ_skew_res[:] = np.NAN

	return IJ_content_res, IJ_skew_res


def init_cumul_nuc(sequence_reference, len_genome, window, dict_nuc):
	"""
	Cumulative of nucleotide count along genome (init from first window)
	"""
	# ATGC
	cumul = np.zeros((4,(len_genome+window) ))

	for j in range(window):
		nuc = sequence_reference[j]
		if nuc in dict_nuc:
			cumul[dict_nuc[nuc]][j] += 1

	return cumul

################################ INPUT DATA ##############################################################################################

# list of files already done
f = open(file_of_filenames_done, 'r')
list_files_done = []
for name in f.readlines():
	list_files_done.append(name.split('\n')[0])
f.close()

print("Already done %d genomes" % len(list_files_done))

# list of all files
f = open(file_of_filenames, 'r')
list_files_genom = []
for name in f.readlines():
	f_name = name.split('\n')[0]
	short_name = f_name.split("/")[-1]
	if short_name not in list_files_done:
		list_files_genom.append(f_name)
f.close()

print("Creating plots for %d genomes" % len(list_files_genom))

#print(list_files_genom)
nb_files = len(list_files_genom)
dict_nuc = {'A' : 0, 'T' : 1, 'G' : 2, 'C' : 3}

for i in range(nb_files):
	# print current file
	filename_ref = list_files_genom[i]
	print("\n%d / %d - %s" %(i, nb_files, filename_ref))

	# read data and init results array
	sequence_reference, len_genome = read_genome(filename_ref)
	if type_window == 'adapt':
		window = int(round(len_genome * frac_window))
	if (len_genome < window) | (window < 1):
		print("Incorrect value for window %d (genome length : %d)" %(window, len_genome))
		continue
	slide_window, seq_right = init_sliding_window(sequence_reference, window)
	GC_content_res, GC_skew_res = init_list_results(len_genome)
	AT_content_res, AT_skew_res = init_list_results(len_genome)
	cumul = init_cumul_nuc(sequence_reference, len_genome, window, dict_nuc)

	# initialisation
	print("Calculating GC skew and AT skew")
	i = 0
	c = collections.Counter(slide_window)
	dict_counts = {'G' : c['G'], 'C' : c['C'], 'A' : c['A'], 'T' : c['T']}
	
	# GC
	sumGC = float(dict_counts['G'] + dict_counts['C'])
	GC_content_res[0][i] = sumGC
	if sumGC > 0:
		GC_skew_res[0][i] = (dict_counts['G'] - dict_counts['C'])/sumGC
	# AT
	sumAT = float(dict_counts['A'] + dict_counts['T'])
	AT_content_res[0][i] = sumAT
	if sumAT > 0:
		AT_skew_res[0][i] = (dict_counts['A'] - dict_counts['T'])/sumAT


	# Compute for all genome
	while(seq_right):
		out_nuc = slide_window.popleft()
		in_nuc = seq_right.popleft()
		slide_window.append(in_nuc)

		i += 1
		if i % 500000 == 0:
			print(i)

		# if in and out are the same : do nothing, append same result
		if out_nuc != in_nuc:
			# remove out from counters
			if out_nuc in dict_nuc:
				dict_counts[out_nuc] -= 1
			if in_nuc in dict_nuc:
				dict_counts[in_nuc] += 1
			sumGC = float(dict_counts['G'] + dict_counts['C'])
			sumAT = float(dict_counts['A'] + dict_counts['T'])

		# fill results
		# GC
		GC_content_res[0][i] = sumGC
		if sumGC > 0:
			GC_skew_res[0][i] = (dict_counts['G'] - dict_counts['C'])/sumGC
		# AT
		AT_content_res[0][i] = sumAT
		if sumAT > 0:
			AT_skew_res[0][i] = (dict_counts['A'] - dict_counts['T'])/sumAT
		# cumul
		if in_nuc in dict_nuc:
			cumul[dict_nuc[in_nuc]][i+window-1] +=1

	GC_content_res = GC_content_res/float(window)
	AT_content_res = AT_content_res/float(window)
	cumul = np.delete(cumul, range(len_genome,cumul.shape[1]),1)
	cumul = np.cumsum(cumul,axis=1)

	Xn = (cumul[dict_nuc['A']] + cumul[dict_nuc['G']]) - (cumul[dict_nuc['C']] + cumul[dict_nuc['T']])
	Yn = (cumul[dict_nuc['A']] + cumul[dict_nuc['C']]) - (cumul[dict_nuc['G']] + cumul[dict_nuc['T']])
	Zn = (cumul[dict_nuc['A']] + cumul[dict_nuc['T']]) - (cumul[dict_nuc['C']] + cumul[dict_nuc['G']])

	GC_content_total = (cumul[dict_nuc['G']][-1] + cumul[dict_nuc['C']][-1]) / float(len_genome)
	AT_content_total = (cumul[dict_nuc['A']][-1] + cumul[dict_nuc['T']][-1]) / float(len_genome)
	ignored_nuc = 1.0 - GC_content_total - AT_content_total

	################################ PLOTS ##############################################################################################
	print("Creating plots")

	pylab.close('all')
	# create figure
	fig, axarr = pylab.subplots(9,1, sharex=True, figsize=(10, 12))

	main_title = "%s - Window size = %d (%.0f %% of genome )\nGC content = %.0f %%, AT content = %.0f %%, ignored = %.0f %%" % (filename_ref.split("/")[-1], window, frac_window*100, GC_content_total*100, AT_content_total*100, ignored_nuc*100)
	pylab.suptitle(main_title, fontsize=16)

	# GC skew
	axarr[0].set_title("GC skew (blue) - Cumulative sum (red)")
	axarr[0].plot(list(GC_skew_res[0]),'b-',alpha=0.5)
	axarr[0].set_ylabel("(G -C) / (G + C)")

	axarr[1].plot(list(np.cumsum(GC_skew_res[0])),'r-',alpha=0.5)
	axarr[1].set_ylabel("(G -C) / (G + C)")
	
	# AT skew
	axarr[2].set_title("AT skew (blue) - Cumulative sum (red)")
	axarr[2].plot(list(AT_skew_res[0]),'b-',alpha=0.5)
	axarr[2].set_ylabel("(A -T) / (A + T)")

	axarr[3].plot(list(np.cumsum(AT_skew_res[0])),'r-',alpha=0.5)
	axarr[3].set_ylabel("(A -T) / (A + T)")

	# Xn
	axarr[4].set_title("Cumulative RY skew (Purine - Pyrimidine)")
	axarr[4].plot(list(Xn),'g-',alpha=0.5)
	axarr[4].set_ylabel("(A + G) - (C + T)")

	# Yn
	axarr[5].set_title("Cumulative MK skew (Amino - Keto)")
	axarr[5].plot(list(Yn),'g-',alpha=0.5)
	axarr[5].set_ylabel("(A + C) - (G + T)")

	# Zn
	axarr[6].set_title("Cumulative H-bond skew (Weak H-bond - Strong H-bond)")
	axarr[6].plot(list(Zn),'g-',alpha=0.5)
	axarr[6].set_ylabel("(A + T) - (G + C)")

	# GC content
	axarr[7].set_title("GC content")
	axarr[7].plot(list(GC_content_res[0]),'k-',alpha=0.5)
	axarr[7].set_ylabel("GC")

	# AT content
	axarr[8].set_title("AT content")
	axarr[8].plot(list(AT_content_res[0]),'k-',alpha=0.5)
	axarr[8].set_ylabel("AT")

	fig.tight_layout()
	fig.subplots_adjust(top=0.88)
	# if save_plots:
	# 	pylab.savefig("%s_%d_%d.png" %(figname_cov,i,i_end))
	# else:
	# 	pylab.show()
	# pylab.show()
	if prefix_figure == 'show':
		pylab.show()
	else :
		pylab.savefig("%s%s_%s_skew.png" %(prefix_figure, filename_ref.split("/")[-1], type_window))
	pylab.clf()
	pylab.close('all')

	###### CLEAR MEMORY
	Xn = None
	Yn = None
	Zn = None
	cumul = None
	# not sure it works


	# # GC skew
	# pylab.figure(figsize=(20, 5))
	# pylab.plot(list(GC_skew_res[0]),'b-',alpha=0.5)
	# pylab.title("GC skew (window size = %d)" %window)
	# pylab.xlabel("Position on reference genome")
	# pylab.ylabel("(G -C) / (G + C)")
	# pylab.show()

	# # AT skew
	# pylab.figure(figsize=(20, 5))
	# pylab.plot(list(AT_skew_res[0]) +list(AT_skew_res[0]),'b-',alpha=0.5)
	# pylab.title("AT skew (window size = %d)" %window)
	# pylab.axvline(x=len_genome, linewidth=1, color='r')
	# pylab.xlabel("Position on reference genome")
	# pylab.ylabel("(A -T) / (A + T)")
	# pylab.show()


"""
for i in range(1,21):
	print("\ninit : %d" % i)
	if i % 5 == 0:
		print("next")
		continue

	if i %11 == 0:
		print(11)
		break
	print("end : %d" %i)

"""