# -*- coding: utf-8 -*-
# arg1 : kraken.out
# arg2 : names_lineage.csv
# arg3 : output.csv

# This takes in input the kraken result (kraken.out) and a csv file with names of species and real taxonomy
# columns of the input csv file must be : "name","full_name","superkingdom","phylum","class","order","family","genus","species"

# this will output a csv file that is a copy of the input, with additionnal columns :
# "good_classification_at_level","wrong_classification_at_level","unknown_taxon_at_level",
# "good_classification_above_level","wrong_classification_above_level","unknown_taxon_above_level","Unclassified"
# each column contains the count of the corresponding classification found by kraken


################################ IMPORT ################################################################################################
import os
import sys
from sequana.lazy import pylab
from sequana.lazy import pandas as pd

scriptpath = "/home/mcardon/Mel/Code/sequana/sequana/"
sys.path.append(os.path.abspath(scriptpath))
import kraken

################################ PARAMETERS ##############################################################################################

filename_kraken = str(sys.argv[1])
filename_read_lineage = str(sys.argv[2])
filename_output = str(sys.argv[3])

level_classification_to_measure = "genus"
# possible choices : ["superkingdom","phylum","class","order","family","genus","species"]

################################ FUNCTIONS ##############################################################################################

def find_tax_info(k_result, tax):
	try:
		# try in records
		info_taxon = k_result.tax.records[tax] 
		#print("Found in record %d" % tax)
	except KeyError:
		try:
			# try in taxonomy
			info_taxon = k_result.tax[tax]
			#print("Found in tax %d" % tax)
		except KeyError:
			#print("Taxon %d not found : fetching from Ensembl" %tax)
			info_taxon = k_result.tax.fetch_by_id(tax)
	return info_taxon


def get_lineage_tax(k_result,info_taxon, res):
	i=0
	if "rank" in info_taxon:
		rank = info_taxon["rank"]
	else:
		rank = "Unknown"

	while rank != "superkingdom":
		try:
			rank = info_taxon["rank"]
			if rank in res:
				res[rank] = info_taxon["scientific_name"]
		except KeyError:
			pass
		
		try:
			tax = int(info_taxon["parent"])
		except TypeError:
			tax = int(info_taxon["parent"]["id"])
		info_taxon  = find_tax_info(k_result, tax)
		#print(info_taxon)
		if i > 500:
			print("Problem in lineage")
			print(info_taxon)
			break
	return res


def classification_result_at_level(df, real_genus, level_classification_to_measure):
	"""
	This calculates the number of good classification, wrong classification, and unknown at this level.
	Output is :
	cl         : list, counts of good, wrong, unknown
	df_unknown : the dataframe of unknown at this level
	"""
	cl = [0]*3
	df_unknown = pd.DataFrame()
	for gen in df[level_classification_to_measure].unique():
		if gen == real_genus:
			#print("identical")
			cl[0] = df[df[level_classification_to_measure] == gen].shape[0]
		elif (str(gen) == str(0)) | (str(gen) == "Unknown"):
			#print("Unknown")
			df_unknown = df[df[level_classification_to_measure] == gen]
			cl[2] = df_unknown.shape[0]
		else:
			#print("wrong")
			cl[1] = df[df[level_classification_to_measure] == gen].shape[0]

	return cl, df_unknown

################################ INPUT DATA ##############################################################################################

k_result = kraken.KrakenResults(filename_kraken)

################################ EXECUTE ##############################################################################################

### get lineage info for all taxons found
df_result = k_result.df
taxons = df_result["taxon"].unique()
result = []
level_tax = ["superkingdom","phylum","class","order","family","genus","species"]

levels_above = []
i = 0
while (level_tax[i] != level_classification_to_measure) & (i < len(level_tax)):
	levels_above.append(level_tax[i])
	i += 1

for tax in taxons:
	# res = [superkingdom,phylum,class,order,family,genus,species]
	res = {"ID":tax, "superkingdom": 0, "phylum": 0, "class": 0, "order": 0, "family": 0, "genus": 0, "species" : 0}
	# try to find taxonomy information
	info_taxon = find_tax_info(k_result, tax)

	# if rank info provided, get lineage
	try:
		#rank = info_taxon["rank"]
		# get lineage here
		#print("Get lineage for taxon %d" % tax)
		res = get_lineage_tax(k_result,info_taxon, res)
	except TypeError:
		res["superkingdom"] = "Unknown"
		#print("Info not found on NCBI for taxon %d" % tax)
	except KeyError:
		pass
		#print("Rank not found for taxon %d" %tax)
		#print(tax, info_taxon.keys())

	result.append(res)


list_result = [ [taxons[i]] + [result[i][lev] for lev in level_tax] for i in range(len(taxons)) ]
df_taxons_lineage = pd.DataFrame(list_result)
df_taxons_lineage.columns = ["taxon"] + level_tax

df_result_merged = df_result.merge(df_taxons_lineage,how="left",on="taxon")

### get species names in kraken.out
names = []
f = open(filename_kraken, 'r')
contig_info = f.readline()
while contig_info:
	contig_list = contig_info.replace("\n","").split("\t")
	name = contig_list[1].replace("_HiSeq","").replace("_MiSeq","").split(".")[0]
	names.append(name)
	contig_info = f.readline()
f.close()
df_result_merged["name"] = names

# get read lineages
df_read_lineage = pd.read_csv(filename_read_lineage)

# check if classification is correct
classification = []
for name in df_read_lineage["name"]:
	df = df_result_merged[ (df_result_merged["status"]=='C') & (df_result_merged["name"] == name)]
	real_genus = df_read_lineage[df_read_lineage["name"] == name].loc[:,level_classification_to_measure].values[0]

	# classified at this level
	cl, df_unknown = classification_result_at_level(df, real_genus, level_classification_to_measure)
	
	# not classified at this level, but classified above
	i = len(levels_above)-1
	lev_ab = levels_above[i]
	cl_above_all = [0]*2
	while (df_unknown.shape[0] > 0) & (i>=0) :
		#print(name, i)
		real_genus = df_read_lineage[df_read_lineage["name"] == name].loc[:,levels_above[i]].values[0]
		cl_above, df_unknown = classification_result_at_level(df_unknown, real_genus, levels_above[i])
		cl_above_all = [ cl_above_all[j] + cl_above[j] for j in range(len(cl_above_all))]
		i -= 1
	cl_above_all.append(cl_above[2])
	cl = cl + cl_above_all
	
	# unclassified
	df = df_result_merged[ (df_result_merged["status"]=='U') & (df_result_merged["name"] == name)]
	cl.append(df.shape[0])
	
	classification.append(cl)

classification_columns = ["good_classification_at_level", "wrong_classification_at_level", "unknown_taxon_at_level",
"good_classification_above_level", "wrong_classification_above_level", "unknown_taxon_above_level","Unclassified"]
classification = pd.DataFrame(classification)
classification.columns = classification_columns

df_read_lineage_result = pd.concat([df_read_lineage,classification],axis=1)
df_read_lineage_result.to_csv(filename_output,index=None)


"""
df_result_merged[(df_result_merged["status"]=='U') & (df_result_merged["name"] == "P_fermentans")]
df_result_merged[(df_result_merged["status"]=='C') & (df_result_merged["name"] == "P_fermentans")]

tax = 610130
info_taxon = find_tax_info(k_result, tax)

res = {"ID":tax, "superkingdom": 0, "phylum": 0, "class": 0, "order": 0, "family": 0, "genus": 0, "species" : 0}
get_lineage_tax(k_result,info_taxon, res)

"""

