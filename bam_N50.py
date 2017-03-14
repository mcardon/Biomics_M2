# -*- coding: utf-8 -*-
# arg1 : bam file
# arg2 : type of filter : 'length' , 'fraction' , 'no'
# arg3 : if length : int (minimum length to keep read), if fraction : float (fraction of longest reads to keep, between 0 and 1) , if no : nothing (no output)
# arg4 : bam output file

# write a bam file filtered

################################ IMPORT ################################################################################################

import sys
import pysam


################################ PARAMETERS ##############################################################################################

filename_BAM = str(sys.argv[1])
type_filter = str(sys.argv[2])

if len(sys.argv) > 3:
	if type_filter == 'length':
		min_len = int(sys.argv[3])
		filename_output = str(sys.argv[4])
	elif type_filter == 'fraction':
		frac_keep = float(sys.argv[3])
		filename_output = str(sys.argv[4])


################################ INPUT DATA ##############################################################################################


b = pysam.AlignmentFile(filename_BAM, check_sq=False)


################################ FUNCTIONS ##############################################################################################


def write_bam_min_len(b, min_len, filename_output):
	"""
	"""
	b_filtered = pysam.AlignmentFile(filename_output, "wb", template=b)
	for read in b:
		if read.query_length >= min_len:
			b_filtered.write(read)
	b_filtered.close()


def set_min_len(b,frac_keep):
	all_read_len = []
	for read in b:
		all_read_len.append(read.query_length)
	# for fraction filtering set min_len
	all_read_len = sorted(all_read_len)
	pos_to_keep = int(round(frac_keep*len(all_read_len)))
	min_len = all_read_len[-pos_to_keep]
	return min_len


################################ EXECUTE ##############################################################################################


# for length filtering : open output
if type_filter == "length":
	# write output
	write_bam_min_len(b, min_len, filename_output)
	b.close()

elif type_filter == "fraction":
	min_len = set_min_len(b,frac_keep)
	b.close()
	# open again
	b = pysam.AlignmentFile(filename_BAM, check_sq=False)
	write_bam_min_len(b, min_len, filename_output)
	b.close()






#b_filtered.close()


