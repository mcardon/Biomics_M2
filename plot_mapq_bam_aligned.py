# -*- coding: utf-8 -*-
# arg1 : BAM file (aligned)
# arg2 : name of output plot
# arg3 : optional threshold for filter (default = 1)

import pysam
import matplotlib.pyplot as plt
import sys
import numpy as np


filename = str(sys.argv[1])
filename_output = str(sys.argv[2])

if len(sys.argv) > 3:
	threshold = int(sys.argv[3])
else:
	threshold = 1


# read input bam
b = pysam.AlignmentFile(filename, mode='rb', check_sq=False)


j = 0
mapq = []

for read in b:
	mapq.append(read.mapq)

	j +=1
	if j%10000 == 0:
		print("Already read %d sequences" %(j))
		#break

plt.hist(mapq, bins=40)
plt.axvline(threshold,color='red')
plt.xlabel("mapq")
plt.yscale('log')
plt.title("reads mapq threshold = %d , N total = %d\n %s" % (threshold , j, filename))
plt.savefig(filename_output)
plt.show()
plt.close()


b.close()

