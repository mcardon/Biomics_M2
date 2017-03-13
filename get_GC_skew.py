# -*- coding: utf-8 -*-
# !! only for single chromosome (bact genome) 


################################ IMPORT ################################################################################################
import pylab
from Bio import SeqIO
import numpy as np
import collections
# import sys

################################ PARAMETERS ##############################################################################################


file_of_filenames.readlines() = str(sys.argv[1])
window = str(sys.argv[2])


# filename_ref = "/home/mcardon/Mel/Datas/Ref_genome_NC/NC_002929_Bpertussis.fasta"
# window = 500000 

################################ FUNCTIONS ##############################################################################################


def init_sliding_window(filename_ref, window):
	"""
	Read fasta file, init sliding window

	return :
	len_genome : int : lenght of reference sequence
	slide_window : deque of n first nucleotides (size of window)
	seq_right : deque of the rest of the genome, with the n-1 nucleotides at the end (circular DNA)
	"""

	# read genome reference
	ref_fasta = SeqIO.parse(filename_ref, "fasta")
	sequence_reference = str(next(ref_fasta).seq).upper()
	len_genome = len(sequence_reference)

	
	# split sequence_reference in two : sliding window and seq_right
	# for circular genome : add begining at the end of genome
	slide_window = collections.deque(sequence_reference[0:window])
	seq_right = collections.deque( sequence_reference[window:len(sequence_reference)] + sequence_reference[0:(window-1)] )

	if len(sequence_reference) < window:
		print("Error : ref genome shorter than window")
		return 0
	else:
		return len_genome, slide_window, seq_right


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



################################ INPUT DATA ##############################################################################################


# list of all files all files
list_files_genom = file_of_filenames.readlines().split('\n')
nb_files = len(list_files_genom)

for i in range(nb_files):
	# print current file
	filename_ref = list_files_genom[i]
	print("%d / %d - %s" %(i, nb_files, filename_ref))

	# read data and init results array
	len_genome, slide_window, seq_right = init_sliding_window(filename_ref, window)
	GC_content_res, GC_skew_res = init_list_results(len_genome)
	AT_content_res, AT_skew_res = init_list_results(len_genome)

	# initialisation
	print("Calculating GC skew")
	i = 0
	c = collections.Counter(slide_window)
	dict_counts = {'G' : c['G'], 'C' : c['C'], 'A' : c['A'], 'T' : c['T']}
	sumGC = float(dict_counts['G'] + dict_counts['C'])

	GC_content_res[0][i] = sumGC
	if sumGC > 0:
		GC_skew_res[0][i] = (G_count - C_count)/sumGC


	while(seq_right):
		out_nuc = slide_window.popleft()
		in_nuc = seq_right.popleft()
		slide_window.append(in_nuc)

		i += 1
		if i % 50000 == 0:
			print(i)

		# if in and out are the same : do nothing
		if out_nuc != in_nuc:
			# remove out from counters
			if out_nuc == 'G':
				G_count -= 1
			elif out_nuc == 'C':
				C_count -= 1

			# add in to counters
			if in_nuc == 'G':
				G_count += 1
			elif in_nuc == 'C':
				C_count += 1

		sumGC = float(G_count + C_count)
		GC_content_res[0][i] = sumGC
		if sumGC > 0:
			GC_skew_res[0][i] = (G_count - C_count)/sumGC


	GC_content_res = GC_content_res/float(window)

	# print(GC_content_res)
	# print(GC_skew_res)

	################################ PLOTS ##############################################################################################
	print("Creating plots")

	pylab.figure(figsize=(20, 5))
	pylab.plot(list(GC_skew_res[0]) +list(GC_skew_res[0]),'b-',alpha=0.5)
	pylab.title("GC skew (window size = %d)" %window)
	pylab.axvline(x=len_genome, linewidth=1, color='r')
	pylab.xlabel("Position on reference genome")
	pylab.ylabel("(G -C) / (G + C)")
	pylab.show()




