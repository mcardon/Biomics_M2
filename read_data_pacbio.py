# -*- coding: utf-8 -*-
# arg1 : BAM file (output of Pacbio sequencer)
import pysam
import matplotlib.pyplot as plt
import sys
import numpy as np
import collections

filename = str(sys.argv[1])

b = pysam.AlignmentFile(filename, check_sq=False)

total_sum = 0
nb_reads = 0

j  = 0

all_GC_percent = []
all_read_len = []
zmw=[]

for read in b:
	# count GC
	c = collections.Counter(read.query_sequence)
	GC_count = c['g'] + c['G'] + c['c'] + c['C']

	rlen = read.query_length
	all_GC_percent.append(GC_count/float(rlen))
	all_read_len.append(rlen)

	# count number of passes for each zmw
	hole_nb = read.qname.split('/')[1]
	zmw.append(hole_nb)
	
	#if break stop before end of file
	j +=1
	if j%10000 == 0:
		print("Already read %d sequences" %(j))
		#print(zmw)
		#break

# close bam file
b.close()



mean_len = np.mean(all_read_len)
mean_GC =  np.mean(all_GC_percent)
print("\nTotal number of reads : %s" % j )
print("Mean read length : %.2f" %(mean_len) )


#all_zmw = zmw.keys()
#distrib_nb_passes = [zmw[h]['nb_passes'] for h in all_zmw]
#distrib_len = [zmw[h]['rlen'] for h in all_zmw]

# distrib_nb_passes = [zmw[h][1] for h in all_zmw]
# distrib_len = [zmw[h][0] for h in all_zmw]

zmw_passes = collections.Counter(zmw)
distrib_nb_passes = [zmw_passes[z] for z in zmw_passes.keys()]

print("Number of passes")
max_nb_pass = max(distrib_nb_passes)
count_pass = collections.Counter(distrib_nb_passes)
for i in range(1,max_nb_pass+1):
	if count_pass[i] != 0:
		print("%d : %d" % (i, count_pass[i]))
print("\n")




######################################## PLOTS ##########################################

# plot GC percent vs read len
plt.plot(all_read_len, all_GC_percent, 'bo', alpha=0.07)
plt.xlabel("Read length")
plt.ylabel("GC percent")
plt.title("GC content vs length \n Mean length = %.2f , Mean GC : %.2f" %(mean_len, mean_GC))
plt.savefig("GC_content_vs_len.png")
#plt.show()
plt.close()

# histogram GC percent
plt.hist(all_GC_percent, bins=40)
plt.xlabel("GC percent")
plt.title("GC content  \n Mean GC : %.2f" %(mean_GC))
plt.savefig("GC_content_hist.png")
#plt.show()
plt.close()

# histogram read length
plt.hist(all_read_len, bins=40)
plt.xlabel("Read length")
plt.title("Read length \n Mean length = %.2f" %(mean_len))
plt.savefig("Read_len_hist.png")
#plt.show()
plt.close()


# histogram nb passes
plt.hist(distrib_nb_passes, bins = max_nb_pass)
plt.xlabel("Nb passes")
plt.yscale('log')
plt.title("Nb passes")
plt.savefig("Nb_passes_hist.png")
#plt.show()
plt.close()






