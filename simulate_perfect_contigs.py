# -*- coding: utf-8 -*-
# arg1 : filename reference genome
# arg2 : filename output.fasta
# arg3 : optionnal : int : begining of the sequence to choose contigs
# arg4 : optionnal : int : length of the sequence (from start chosen) to choose contigs

### draft script : change quotes to execute the part you want

################################ IMPORT ################################################################################################


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sequana.lazy import numpy as np
import sys
import random

################################ PARAMETERS ################################################################################################

filename_ref = str(sys.argv[1])
filename_output = str(sys.argv[2])

if len(sys.argv) > 3:
    start = int(sys.argv[3])
    if len(sys.argv) > 4:
        limit = int(sys.argv[4])
    else:
        print("No limit defined : use the length of the sequence (from start point)")
else:
    start = 0


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

ref = SeqIO.read(filename_ref, "fasta")

"""

# simulate contigs
random_contigs=[]

# if limit not already defined
if not (len(sys.argv) > 4):
    limit=len(ref.seq) - start

while limit > 10000:
    len_contig = int(round(random.gauss(100000,100000)))
    if (len_contig > 0) & (len_contig < limit):
        end = start+len_contig
        contig = ref.seq[start:end]
        record = SeqRecord(contig,'contig_%i_%i_len_%i' % (start, end,end-start), '', '')
        random_contigs.append(record)
        limit = limit +(start -end)
        start = end


SeqIO.write(random_contigs, filename_output, "fasta")
"""

"""
#### create fasta with 2 times the reference sequence
sequence_reference, len_genome = read_genome(filename_ref)
record_2x = SeqRecord(Seq(sequence_reference + sequence_reference,generic_dna), id = "NC_002929 genome duplicated")
SeqIO.write(record_2x, "NC_002929_duplicated.fasta", "fasta")
"""


# create N rotations of reference genome
sequence_reference, len_genome = read_genome(filename_ref)
seq_duplicated = sequence_reference+sequence_reference

N = 10
start_point = np.linspace(0,len_genome,N+1)
start_point = [int(round(i)) for i in start_point[0:(len(start_point)-1)]]

leak = start_point[1] # number of bases to keep at the end of the rotation : duplicate of the begining
print(leak)

for start in start_point:
    seq_rotated = seq_duplicated[start:(start+len_genome+leak)]
    record_rotated = SeqRecord(Seq(seq_rotated,generic_dna), id = "NC_002929 genome rotated start=" +str(start))
    SeqIO.write(record_rotated, "NC_002929_rotate_"+ str(start) +".fasta", "fasta")



#len_genome = 4086189






