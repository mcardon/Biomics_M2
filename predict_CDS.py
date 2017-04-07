# -*- coding: utf-8 -*-
# !! only for single chromosome (bact genome) fasta must contain only one entry
# arg1 : filename genome (fasta format)
# arg2 : prefix filename for plots and csv file (if no : don't save anything, show plots)
# arg3 : threshold : int : minimum size to output the ORF or CDS (DNA sequence size)
# arg4 : type of sequence to witch the threshold is applied : possible choices : ORF CDS

# Input sequence is assumed to be 5'->3' and with only ATGC characters
# Output position of ORF and CDS is 0-based

################################ IMPORT ################################################################################################
import pylab
from Bio import SeqIO
import numpy as np
import pandas as pd
import sys
import time


################################ PARAMETERS ##############################################################################################


filename_fasta = str(sys.argv[1])
filename_output = str(sys.argv[2])
THRESHOLD = int(sys.argv[3])
TYPE_FILTER = str(sys.argv[4])

if TYPE_FILTER not in ["ORF","CDS"]:
	print("Invalid arg4 : type filter. Possible choices : ORF CDS\nConsidering CDS size by default")

save_output = True
if filename_output == "no":
	save_output = False

CODONS_STOP =     ["TAA","TGA","TAG"]
CODONS_STOP_REV = ["TTA","TCA","CTA"]
CODONS_START =     ["ATG"]
CODONS_START_REV = ["CAT"]

################################ FUNCTIONS ##############################################################################################

def read_genome(filename_fasta):
    """
    """
    # read genome reference
    ref_fasta = SeqIO.parse(filename_fasta, "fasta")
    sequence_reference = str(next(ref_fasta).seq).upper()
    len_genome = len(sequence_reference)
    return sequence_reference, len_genome



################################ FUNCTIONS ##############################################################################################

def update_ORF_frame(i, frame, codon, ORF_pos, begin_f, end_f, begin_r, end_r, pos_ATG_f, pos_ATG_r):
		codon = codon + nuc

		if len(codon) == 3:
			# test for start forward (if none already found)
			is_start = codon in CODONS_START
			if is_start & np.isnan(pos_ATG_f):
				pos_ATG_f = i-2
			# test for stop forward
			is_stop = codon in CODONS_STOP
			if is_stop:
				end_f = i
				# test is length of CDS or ORF found is enough to be in output
				if TYPE_FILTER == "CDS":
					len_to_filter = end_f - pos_ATG_f # len_CDS
				else:
					len_to_filter = end_f - begin_f   # len_ORF

				if len_to_filter > THRESHOLD:
					len_ORF = end_f - begin_f
					len_CDS = end_f - pos_ATG_f
					ORF_pos.append([begin_f,end_f,frame+1, pos_ATG_f, len_ORF, len_CDS])
				begin_f = i+1
				pos_ATG_f = np.nan

			# test for start reverse
			is_start_rev = codon in CODONS_START_REV
			if is_start_rev :
				pos_ATG_r = i
			# test for stop reverse
			is_stop_rev = codon in CODONS_STOP_REV
			if is_stop_rev:
				end_r = i-3
				# test is length of CDS or ORF found is enough to be in output
				if TYPE_FILTER == "CDS":
					len_to_filter = pos_ATG_r - begin_r # len_CDS
				else:
					len_to_filter = end_r - begin_r   # len_ORF

				if len_to_filter > THRESHOLD:
					len_ORF = end_r - begin_r
					len_CDS = pos_ATG_r - begin_r
					ORF_pos.append([begin_r,end_r,-(frame+1), pos_ATG_r, len_ORF, len_CDS])
				begin_r = i-3+1
				pos_ATG_r = np.nan

			# reset codon
			codon = ""

		return codon, ORF_pos, begin_f, end_f, begin_r, end_r, pos_ATG_f, pos_ATG_r




################################ EXECUTE ##############################################################################################


#sequence_reference = "atggattaaccccgcgaatgatgtgatcgtcgttggtggcggtcacgccggtacggaggcagccctggctgcagcccgcgccggcgcacagacattgctgcttacccacaatatcgag".upper()

time_begin = time.time()
sequence_reference, len_genome = read_genome(filename_fasta)

print("Sequence size : %d" %len_genome)

# init variables
codon = ["","-","--"]
begin_f = [0]*3
begin_r = [0]*3
end_f = [0]*3
end_r = [0]*3
pos_ATG_f = [np.nan]*3
pos_ATG_r = [np.nan]*3

# result list
ORF_pos = []

# search ORF and CDS
for i, nuc in enumerate(sequence_reference):
	#print(i, nuc)
	frame = (i+3-2) % 3
	if (i % 200000)==0:
		print("%d / %d" %(i, len_genome))
	for j in range(3):
		codon[j], ORF_pos, begin_f[j], end_f[j], begin_r[j], end_r[j], pos_ATG_f[j], pos_ATG_r[j] = update_ORF_frame(i, frame, codon[j], ORF_pos, begin_f[j], end_f[j], begin_r[j], end_r[j], pos_ATG_f[j], pos_ATG_r[j])


# convert result to dataframe and output csv file
df = pd.DataFrame(ORF_pos)
df.columns = ["begin_pos","end_pos","frame","pos_ATG","len_ORF","len_CDS"]
filename_output_csv = filename_output + "_" + TYPE_FILTER + "_" + str(THRESHOLD) + ".csv"
if save_output:
	df.to_csv(filename_output_csv, index=False)

time_end = time.time()


print("Finding %s : Finished in %.2f seconds, plotting result" % (TYPE_FILTER,(time_end - time_begin)))

n_ORF = df["len_ORF"].dropna().shape[0]
n_CDS = df["len_CDS"].dropna().shape[0]

# plot for all ORF and CDS
pylab.hist(df["len_ORF"].dropna(),alpha=0.5, label="ORF, N = " + str(n_ORF),bins=40)
pylab.hist(df["len_CDS"].dropna(),alpha=0.5, label="CDS, N = " + str(n_CDS),bins=40)
pylab.xlabel("Length")
pylab.ylabel("#")
pylab.legend()
pylab.title("Length of ORF and CDS (after filter %s > %d)" %(TYPE_FILTER, THRESHOLD))
if save_output:
	pylab.savefig(filename_output + "_hist.png")
else:
	pylab.show()


pylab.hist(df["len_ORF"].dropna(),alpha=0.5, label="ORF, N = " + str(n_ORF),bins=40)
pylab.hist(df["len_CDS"].dropna(),alpha=0.5, label="CDS, N = " + str(n_CDS),bins=40)
pylab.xlabel("Length")
pylab.ylabel("#")
pylab.yscale("log")
pylab.legend()
pylab.title("Length of ORF and CDS (after filter %s > %d)" %(TYPE_FILTER, THRESHOLD))
if save_output:
	pylab.savefig(filename_output + "_hist_logscale.png")
else:
	pylab.show()

# number of ORF and CDS found by frame
frames = [-3, -2, -1, 1, 2, 3]
nb_res_ORF = []
nb_res_CDS = []
for fr in frames:
	nb_res_ORF.append(df[df["frame"] == fr]["len_ORF"].dropna().shape[0])
	nb_res_CDS.append(df[df["frame"] == fr]["len_CDS"].dropna().shape[0])

bar_width = 0.35
pylab.bar(np.array(frames)-(bar_width/2), nb_res_ORF, bar_width, alpha=0.5, label="ORF")
pylab.bar(np.array(frames)+(bar_width/2), nb_res_CDS, bar_width, alpha=0.5, label="CDS")
pylab.legend(loc=1)
pylab.title("Number of ORF and CDS by frame")
pylab.show()


print("End")


"""
for orf in ORF_pos:
	print(orf)
	print(sequence_reference)
	print("-"*orf[0] + sequence_reference[orf[0]:(orf[1]+1)])
	if np.isnan(orf[3]):
		print("No start")
	else:
		print(" "*orf[3] + "*")

"""




