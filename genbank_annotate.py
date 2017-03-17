# -*- coding: utf-8 -*-
# arg1 : input.csv : file to annotate in csv format
# arg2 : position_col_name : str : name of the column containing position information
# arg3 : output.csv
# arg4 : genbank_file.gbk


# this will add genbank annotation in the last columns and output in csv format
# this will parse only the first entry in genbank file (first genome) so better use single genome genbank file


##### Bug : if key not found in features, error
##### check if key inside dict before


################################ IMPORT ################################################################################################
import sys
from sequana.lazy import pandas as pd
from Bio import SeqIO

################################ PARAMETERS ################################################################################################

file_input = str(sys.argv[1])
col_name_pos = str(sys.argv[2])
file_output = str(sys.argv[3])
file_genbank = str(sys.argv[4])

################################ INPUT DATA ################################################################################################

df = pd.read_csv(file_input,sep=",")
df = df.sort_values(by=col_name_pos,ascending=True)

seq_records = SeqIO.parse(file_genbank, "genbank")
record = next(seq_records)


# for each position, check if there is CDS annotation
rec_i = 0
b = 0
e = 0

# init table all results
result_annot = []
header_df_results = ["CDS"," start"," end"," strand"," gene_ID"," gene_name"," product"," note"]

for pos in list(df['position']):

	# find CDS just after variant
	while(pos > e):
		rec_i += 1
		# if CDS found : update position of last CDS seen
		if record.features[rec_i].type == "CDS":
			b = record.features[rec_i].location.nofuzzy_start
			e = record.features[rec_i].location.nofuzzy_end
	
	if (pos >= b) & (pos <= e):
		########### variant inside CDS
		res = []
		# convert strand info
		strand = '+' if record.features[rec_i].strand > 0 else '-'
		# CDS, start, end, strand, gene_ID, gene_name, product, note
		res.append(record.features[rec_i].type)                               # CDS
		res.append(b)                                                         # start
		res.append(e)                                                         # end
		res.append(strand)                                                    # strand
		quals = record.features[rec_i].qualifiers
		res.append(" ; ".join(quals["db_xref"]) if "db_xref" in quals else None ) # gene ID (maybe more than one)
		res.append(" ; ".join(quals["gene"]) if "gene" in quals else None )       # gene name (maybe more than one)
		res.append(" ; ".join(quals["product"]) if "product" in quals else None ) # product of the gene
		res.append(" ; ".join(quals["note"]) if "note" in quals else None )       # note (description)
	else:
		########### variant not in CDS : append empty line
		res = [None]*len(header_df_results)

	# append to final result
	result_annot.append(res)

df_result_annot = pd.DataFrame(result_annot)
df_result_annot.columns = header_df_results
df_result_annot.index = df.index

df_ouput = pd.concat([df, df_result_annot], axis = 1)
df_ouput.to_csv(file_output)


################################  ################################################################################################
################################  ################################################################################################
################################  ################################################################################################




