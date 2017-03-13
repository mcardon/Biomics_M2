# -*- coding: utf-8 -*-
################################ IMPORT ################################################################################################
import collections
import numpy as np
import pylab

################################ PARAMETERS ##############################################################################################

## in the pileup file the position is 1-base
## keep this base in the np array in the loop
## after the loop > convert to 0-base
file_pileup = "/home/mcardon/Mel/Datas/Toy_data_pacbio/toy_data_mapped_blasr_9X_mpileup.pileup"
len_genome = 4100000 # size of genome or larger (unfilled positions will be ignored)

################################ INPUT DATA ##############################################################################################



################################ ACCURACY ##############################################################################################

accuracy_no_ins = np.empty((1,len_genome))
accuracy_no_ins[:] = np.NAN

accuracy_all = np.empty((1,len_genome))
accuracy_all[:] = np.NAN

print("Parsing pileup file")
with open(file_pileup, 'r') as f:
    for line in f:
        l = line.replace("\n","").split("\t")
        if len(l) > 4:
            pos = int(l[1])
            ref = l[2]
            cov = int(l[3])
            bas = l[4]

            # count matches
            bas_count = collections.Counter(bas)
            match = bas_count[','] + bas_count['.']
            # count deletion
            deletion = bas_count['*']

            # delete match and deletion (already counted)
            bas_save = bas
            bas = bas.replace(',','').replace('.','').replace('*','')

            # count insertion and mismatch
            indel = False
            stride = []
            is_insert = False
            count_ins = 0
            mismatch = 0
            if len(bas) > 0:
                i = 0
                while i < len(bas):
                    if indel & bas[i].isdigit():
                        stride.append(bas[i])
                    elif indel & (not bas[i].isdigit()):
                        stride = int(''.join(stride))
                        i += stride-1
                        indel = False
                        if is_insert:
                            count_ins += stride
                            is_insert = False
                        
                    else:
                        if bas[i] == '+':
                            indel = True
                            stride = []
                            is_insert = True
                            
                        elif bas[i] == '-':
                            indel = True
                            stride = []

                        else:
                            mismatch +=1
                    i += 1

                #print("matches = %d, mismatch = %d, del = %d, ins = %d" %(match, mismatch, deletion, count_ins))

                acc = match / float(match + mismatch + deletion)
                accuracy_no_ins[0][pos] = acc

                acc = match / float(match + mismatch + deletion + count_ins)
                accuracy_all[0][pos] = acc

################################ PLOTS ##############################################################################################
print("Creating plots")

pylab.figure(figsize=(20, 5))
pylab.plot(list(accuracy_no_ins[0]),'b.',alpha=0.1)
pylab.title("Accuracy (without insertions)")
pylab.xlabel("Position on reference genome")
pylab.ylabel("Accuracy (without insersion)")
pylab.show()
pylab.clf()

pylab.figure(figsize=(20, 5))
pylab.plot(list(accuracy_all[0]),'b.',alpha=0.1)
pylab.title("Accuracy")
pylab.xlabel("Position on reference genome")
pylab.ylabel("Accuracy")
pylab.show()
















