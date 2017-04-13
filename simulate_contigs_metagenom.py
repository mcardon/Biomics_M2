# -*- coding: utf-8 -*-
# arg1 : fof database of reference genomes (file of filenames)
# arg2 : filename output.fastq
# arg3 : mean contig length to simulate
# arg4 : standard deviation of contig length to simulate
# arg5 : number of contigs to simulate from each reference

## this generates contigs form reference sequences and output as fastq format
## warning : no quality values ! all set to "!", fastq format is only for compatibility with other tools 

################################ IMPORT ################################################################################################

from Bio import SeqIO
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sequana.lazy import numpy as np
import sys
import os
import random

################################ PARAMETERS ################################################################################################

fof_database    = str(sys.argv[1])
filename_output = str(sys.argv[2])
len_mu          = int(sys.argv[3])
len_sigma       = int(sys.argv[4])
nb_contigs      = int(sys.argv[5])

prefix_tmp_files = "./temp_dir_contig_sim/"
suffix_tmp_files = "_simulated_contigs_tmp.fasta"


################################ FUNCTIONS ################################################################################################

def list_files_fof(fof_database):
    """list of all files in fof"""
    f = open(fof_database, 'r')
    list_files = []
    for name in f.readlines():
        f_name = name.split('\n')[0]
        list_files.append(f_name)
    f.close()
    return list_files


def read_genome(filename_ref):
    """read genome reference"""
    ref_fasta = SeqIO.read(filename_ref, "fasta")
    sequence_reference = str(ref_fasta.seq).upper()
    len_genome = len(sequence_reference)
    return sequence_reference, len_genome


################################ EXECUTE ################################################################################################


list_files = list_files_fof(fof_database)

output_handle = open(filename_output, "w")

for filename_ref in list_files:
    sequence_reference, len_genome = read_genome(filename_ref)
    filename_ref_short = filename_ref.split("/")[-1]

    # randoms length and start point for contigs
    len_contigs = [int(abs(round(random.gauss(len_mu,len_sigma))) % len_genome) for i in range(nb_contigs)]
    start_contigs = random.sample(range(len_genome), nb_contigs)
    random_contigs = []

    for i in range(nb_contigs):
        start = start_contigs[i]
        end = start + len_contigs[i]
        #print(start, end)
        if end <= len_genome:
            #print(sequence_reference[start:end])
            contig = sequence_reference[start:end]
        else:
            #print(sequence_reference[start:len_genome])
            #print(sequence_reference[0:(end % len_genome)])
            contig = sequence_reference[start:len_genome] + sequence_reference[0:(end % len_genome)]

        id_contig = "%s_start_%i_len_%i" % (filename_ref_short, start, len_contigs[i])
        qual_contig = "!"*len_contigs[i]
        fastq_string = "@%s\n%s\n+\n%s\n" % (id_contig, contig, qual_contig)

        #record = SeqRecord(Seq(contig, generic_dna ), ,filename_ref,'')
        record = SeqIO.read(StringIO(fastq_string), "fastq")
        random_contigs.append(record)

    SeqIO.write(random_contigs, output_handle , "fastq")

output_handle.close()



