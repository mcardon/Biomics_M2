# -*- coding: utf-8 -*-
# arg1 : BAM file (after alignment)
# arg2 : str : type of alignement (BWA or BLASR)
# arg3 : int : minimum mapq (optionnal, default 0)


####################################################################################
#                                  IMPORT
####################################################################################
import pysam
import matplotlib.pyplot as plt
import sys
import numpy as np

####################################################################################
#                                    ARGS
####################################################################################
filename = str(sys.argv[1])
type_align = str(sys.argv[2])

if len(sys.argv) > 3:
	threshold = int(sys.argv[3])
else:
	threshold = 0


####################################################################################
#                                FUNCTIONS
####################################################################################

def compute_concordance(read, type_align):
	"""
	read : pysam mapped read
	type_align : str : BWA or BlASR

	returns concordance between mapped read and ref
	from cigar data in pysam read
	formula : 1 - (in + del + mismatch / (in + del + mismatch + match) )
	"""

	# concordance from cigar stats
	# [M, I, D, N, S, H, P, =, X, NM]
	# concordance = 1 - (sum(I,D,X) / sum(I,D,X,=))

	if type_align == 'BLASR':
		diff_cig = sum( [read.get_cigar_stats()[0][i] for i in [1,2,8]] )
		id_cig = read.get_cigar_stats()[0][7]
		cig =  1 - (diff_cig / float(diff_cig + id_cig))

	elif type_align == 'BWA':
		cig_stats = read.get_cigar_stats()[0]
		indel_cig = sum( [cig_stats[i] for i in [1,2]] )
		mismatch = cig_stats[10] - indel_cig
		diff_cig = indel_cig + mismatch
		tot = cig_stats[0] + indel_cig
		cig =  1 - (diff_cig / float(tot))

	else:
		print("Incorrect read type. Possible type : BLASR or BWA")
		sys.exit(1)

	return cig



####################################################################################
#                                  EXECUTION
####################################################################################

b = pysam.AlignmentFile(filename, check_sq=False)

concordance = []
mapq = []
mapq_filter = []
read_len = []
j = 0

for read in b:
	# mapping quality
	mapq.append(read.mapq)

	if (not read.is_unmapped) & (read.mapq >= threshold):
		# concordance
		cig = compute_concordance(read, type_align)
		concordance.append(cig)

		# read length
		read_len.append(read.rlen)

		# mapq
		mapq_filter.append(read.mapq)

		#if break stop before end of file
		j +=1
		if j%10000 == 0:
			print("Already read %d sequences" %(j))
			#break

# close bam file
b.close()

# sort list mapq and concordance
#mapq_conc = sorted(zip(mapq_filter,concordance), key=lambda pair: pair[0])
#mapq_filter = [m for (m,c) in mapq_conc]
#concordance = [c for (m,c) in mapq_conc]
levels_mapq = sorted(set(mapq_filter))

data_boxplot_mapq_conc = [ [concordance[i] for i in range(len(mapq_filter)) if mapq_filter[i] == m ] for m in levels_mapq]

####################################################################################
#                                  PRINT RESULTS
####################################################################################

# mean and median
conc_mean = np.mean(concordance)
conc_med = np.median(concordance)
len_mean = np.mean(read_len)
len_med = np.median(read_len)
print("\nTotal number of reads after quality filter : %s" % j )
print("Mean concordance   : %.2f" %(conc_mean) )
print("Median concordance : %.2f" %(conc_med) )
print("Mean length        : %.2f" %len_mean)
print("Median length      : %.2f" %len_mean)



####################################################################################
#                                   PLOTS
#####################################################################################

# histogram concordance
plt.hist(concordance, bins = 40)
plt.xlabel("Concordance")
#plt.yscale('log')
plt.title("Concordance (N = %d)\nMapq > %d , Mean = %.2f , Median = %.2f" %(j,threshold,conc_mean,conc_med))
plt.savefig("Concordance_hist.png")
plt.show()
plt.close()

# histogram mapping quality
plt.hist(mapq, bins=20)
plt.axvline(threshold,color='red')
plt.xlabel("mapq")
plt.yscale('log')
plt.title("reads mapq - threshold = %d , N after filter = %d\n %s" % (threshold , j, filename))
plt.savefig("Mapq.png")
plt.show()
plt.close()

# concordance vs read length
plt.plot(read_len, concordance, 'bo', alpha=0.07)
plt.xlabel("Read length")
plt.title("Concordance vs Read length - Mapq > %d , N after filter = %d\n %s" % (threshold , j, filename))
plt.savefig("Concordance_vs_read_len.png")
plt.show()
plt.close()

# concordance vs mapping quality
plt.figure(figsize=(20, 5))
plt.boxplot(data_boxplot_mapq_conc)
plt.xticks(range(1,len(levels_mapq)+1,1),levels_mapq)
plt.xlabel("Mapq")
plt.ylabel("Concordance")
plt.title("Concordance vs Mapq - Mapq > %d , N after filter = %d\n %s" % (threshold , j, filename))
plt.savefig("Concordance_vs_mapq.png")
plt.show()
plt.close()
