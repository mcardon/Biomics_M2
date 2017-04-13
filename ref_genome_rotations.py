# -*- coding: utf-8 -*-
# arg1 : filename reference genome
# arg2 : filename prefix for rotations output
# arg3 : number of rotations

### fasta input file must contain only one sequence

################################ IMPORT ################################################################################################


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sequana.lazy import numpy as np
import sys
#import random

################################ PARAMETERS ################################################################################################

filename_ref = str(sys.argv[1])
filename_output = str(sys.argv[2])
N = int(sys.argv[3])

################################ FUNCTIONS ################################################################################################

def read_genome(filename_ref):
    """
    """
    # read genome reference
    ref_fasta = SeqIO.read(filename_ref, "fasta")
    sequence_reference = str(ref_fasta.seq).upper()
    len_genome = len(sequence_reference)
    return sequence_reference, len_genome



################################ EXECUTE ################################################################################################


# create N rotations of reference genome
sequence_reference, len_genome = read_genome(filename_ref)
seq_duplicated = sequence_reference+sequence_reference

start_point = np.linspace(0,len_genome,N+1)
start_point = [int(round(i)) for i in start_point[0:(len(start_point)-1)]]

leak = start_point[1] # number of bases to keep at the end of the rotation : duplicate of the begining
#print(leak)

for start in start_point:
    seq_rotated = seq_duplicated[start:(start+len_genome+leak)]
    record_rotated = SeqRecord(Seq(seq_rotated,generic_dna), id = "%s genome rotated start=%s" % (filename_output.split("/")[-1],str(start)))
    SeqIO.write(record_rotated, "%s_rotate_%s.fasta" % (filename_output,str(start)), "fasta")




