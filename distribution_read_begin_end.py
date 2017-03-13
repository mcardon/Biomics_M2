import pysam
import pylab
import matplotlib.pyplot as plt
import numpy as np
import collections



b = pysam.AlignmentFile("m54091_161109_200101.subreads_ref_NC_002929_blasr.bam", check_sq=False)


start = []
end = []
N = 0
for read in b:
    start.append(read.reference_start)
    end.append(read.reference_end)
    N +=1
    if N%10000 == 0:
        print("Already read %d sequences" %(N))
        #break

print("Counting start positions")
c_start = collections.Counter(start)
print("Counting end positions")
c_end =  collections.Counter(end)

print("Find max")
max_pos = max(max(c_start.keys()), max(c_end.keys()))
print(max_pos)
print("Compute values")
k = range(1,max_pos+1)
val_start = [c_start[i] for i in k]
val_end = [c_end[i] for i in k]

# params
fontsize = 12
grid = True

print("Histograms")
# histogram start positions
pylab.plot(k,val_start, 'bo', alpha=0.05)
pylab.xlabel("Position first mapped", fontsize=fontsize)
pylab.ylabel("#", fontsize=fontsize)
pylab.title("Number of reads mapped at starting position",fontsize=fontsize)
if grid is True:
    pylab.grid(True)
pylab.savefig("Mapping_start_pos.png")
pylab.show()

# histogram start positions
pylab.clf()
pylab.plot(k, val_end,'bo', alpha=0.05)
pylab.xlabel("Position last mapped", fontsize=fontsize)
pylab.ylabel("#", fontsize=fontsize)
pylab.title("Number of reads mapped at ending position",fontsize=fontsize)
if grid is True:
    pylab.grid(True)
pylab.savefig("Mapping_end_pos.png")
pylab.show()
pylab.close()