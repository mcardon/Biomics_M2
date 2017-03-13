# -*- coding: utf-8 -*-
# arg1 : reference fasta file
# arg2 : bam aligned file sorted (will raise error if not sorted)

# !! only for single chromosome (bact genome) 



################################ IMPORT ################################################################################################
import pylab
import pysam
from Bio import SeqIO
import collections
import sys

################################ PARAMETERS ##############################################################################################

## input file must be sorted !
## samtools sort -o output_sorted.bam input_unsorted.bam

# filename_ref = str(sys.argv[1])
# filename = str(sys.argv[2])

filename_ref = "/home/mcardon/Mel/Datas/Ref_genome_NC/NC_002929_Bpertussis.fasta"
filename = "/home/mcardon/Mel/Datas/Toy_data_pacbio/toy_data_mapped_blasr_9X_sorted.bam"



################################ INPUT DATA ##############################################################################################

b = pysam.AlignmentFile(filename, check_sq=False)

# read genome reference
ref_fasta = SeqIO.parse(filename_ref, "fasta")
sequence_reference = str(next(ref_fasta).seq).upper()

if not b.has_index():
	pysam.index(filename)
	b = pysam.AlignmentFile(filename, check_sq=False)


################################ ACCURACY ##############################################################################################



# from easydev import do_profile

# @do_profile()
def get_acc():
	i = 0
	covered_positions = []
	all_accuracy      = []
	for p in b.pileup():
		# position in genome
		covered_positions.append(p.pos)

		# ref sequence
		base_ref = sequence_reference[p.pos]
		# get all reads mapped at this position
		bases = []
		for pileupread in p.pileups:
			# get information from read
			if pileupread.is_del:
				bases.append('del')
			# elif pileupread.is_refskip:
			# 	# XXX how to deal with insert in accuracy ?
			else:
				q_pos = pileupread.query_position
				nuc = pileupread.alignment.query_sequence[q_pos]
				bases.append(nuc)
		
		count_bases = collections.Counter(bases)

		accuracy = count_bases[base_ref] / float(sum(count_bases.values()))

		all_accuracy.append(accuracy)

		i +=1
		if i % 10000 == 0:
			print("Read %d positions" %i)
			break
	return covered_positions, all_accuracy


################################ PLOTS ##############################################################################################
covered_positions, all_accuracy = get_acc()
pylab.figure(figsize=(20, 5))
pylab.plot(covered_positions, all_accuracy,'b.',alpha=0.1)
pylab.title("Accuracy")
pylab.xlabel("Position on reference genome")
pylab.ylabel("Accuracy")
pylab.show()


"""
i = 0
for p in b.pileup():
	i +=1
	# position in genome
	print("\ncoverage at base %s = %s" % (p.pos, p.n))
	# ref sequence
	print("Reference = %s" %sequence_reference[p.pos])
	# get all reads mapped at this position
	ref_nuc = []
	for pileupread in p.pileups:
		
		# get reference data (and check it's always the same)
		# it does not work ! (not always the same !!)
		# https://github.com/pysam-developers/pysam/issues/225
		#for alignment_tuple in pileupread.alignment.get_aligned_pairs(with_seq=True):
		#	if alignment_tuple[1] == p.pos:
		#		ref_nuc.append(alignment_tuple[2])
		


		# get information from read
		if pileupread.is_del:
			print("base in read %s = del" %(pileupread.alignment.query_name))
		elif pileupread.is_refskip:
			print("base in read %s = ins"% (pileupread.alignment.query_name))
		else:
			print("base in read %s = %s " %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position] ) )

	if i > 100:
		break

	ref_nuc_set = set(ref_nuc)
	if len(ref_nuc_set) > 1:
		print("Error in getting reference")
		print(ref_nuc)

"""