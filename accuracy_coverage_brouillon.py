# -*- coding: utf-8 -*-
# arg1 : reference fasta file
# arg2 : list.txt : list of filenames of csv coverage files
# arg3 : int : window size
# !! only for single chromosome (bact genome) 


#import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
#import collections
#from sequana import *


filename_ref = str(sys.argv[1])
list_filenames_coverage = str(sys.argv[2])
window = int(sys.argv[3])


data = pd.DataFrame()

# read genome reference
ref_fasta = SeqIO.parse(filename_ref, "fasta")

f = open(list_filenames_coverage, 'r')
l = f.readline()

# read all files
while(l):
	file_cov = l.split('\n')[0]
	data_cov = pd.read_csv(file_cov, sep=',', header=None)
	data = pd.concat([data, data_cov] ,axis=0)
	print(data.shape)
	l = f.readline()
f.close()


# sort datahe
data.sort_values(by = 0, axis=0, inplace=True)
data.columns = ['pos','pos1','coverage']
# Coverage data : 1 based coordinates
# Fasta ref     : 0 based corrdinates
# convert coverage data into 0 based coord
data.index = data['pos'] -1
data.drop('pos', axis=1,inplace=True)
data.drop('pos1', axis=1,inplace=True)

# add col for GC content
data['ins'] = np.nan
data['del'] = np.nan
data['mismatch'] = np.nan

print(data.iloc[0:5,:])
print(data.iloc[0,:])







