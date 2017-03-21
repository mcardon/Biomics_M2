# -*- coding: utf-8 -*-
######### Do not work (because it would take too long time, useless to finish this)
# !! only for single chromosome (bact genome) 
# arg1 : contigs.fasta
# arg2 : reference.fasta

# works with skbio 0.5.1

################################ IMPORT ################################################################################################
from sequana.lazy import pylab
from sequana.lazy import pandas as pd
from skbio.alignment import StripedSmithWaterman
from Bio import SeqIO
import sys
import time

print("Import done")

################################ PARAMETERS ################################################################################################

filename_contigs = str(sys.argv[1])
filename_ref = str(sys.argv[2])




################################ FUNCTIONS ################################################################################################


def read_genome(filename_ref):
	"""
	"""
	# read genome reference
	ref_fasta = SeqIO.parse(filename_ref, "fasta")
	sequence_reference = str(next(ref_fasta).seq).upper()
	len_genome = len(sequence_reference)
	return sequence_reference, len_genome


################################ EXECUTE ################################################################################################

# short test
# query = StripedSmithWaterman("ACTAAGGCTCTCTACCCCTCTCAGAGA")
# alignment = query("AAAAAACTCTCTAAACTCACTAAGGCTCTCTACCCCTCTTCAGAGAAGTCGA")
# print(alignment)
sequence_reference, len_genome = read_genome(filename_ref)
sequence_reference_2x = StripedSmithWaterman(sequence_reference + sequence_reference)
print("reference sequence loaded")

x = range(100,10000,300)
res_time = []
for i in x:
	time1 = time.clock()
	alignment_best = sequence_reference_2x(sequence_reference[0:i])
	time2 = time.clock()
	print("Align %d pb in %f, score = %d" % (i,time2 - time1,alignment_best.optimal_alignment_score))
	res_time.append(time2 - time1)

df_res_time = pd.DataFrame([x,res_time])
df_res_time = df_res_time.transpose()
df_res_time.columns = ["len_seq", "time_minutes"]
df_res_time["time_minutes"] = df_res_time["time_minutes"] / float(60)
df_res_time.to_csv("2017_03_20_SW_time.csv")



pylab.plot(df_res_time["len_seq"], df_res_time["time_minutes"],"b-")
pylab.xlabel("Length of aligned sequence")
pylab.ylabel("Time (minutes)")
pylab.title("Time for Smith waterman alignment with scikit 0.5.1")
pylab.show()

estim = df_res_time["time_minutes"].iloc[-1]*len_genome / float(df_res_time["len_seq"].iloc[-1]*60) # in minutes

print("Time estimation (linear) for contig of %d = %f minutes (%f hours)" %(len_genome, estim, estim/60))





