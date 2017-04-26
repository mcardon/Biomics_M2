# -*- coding: utf-8 -*-
# arg1 : kraken.out
# arg2 : output.csv

## Plot section to finish (or do it in another script)

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
filename_output = str(sys.argv[2])


################################ FUNCTIONS ##############################################################################################

def parse_details_kraken(detail_kr):
	detail_dict = {}
	for kmer_found in detail_kr:
		kmer_found = kmer_found.split(":")
		try:
			kmer_found = [int(k) for k in kmer_found]
		except ValueError:
			print(detail_kr)
			print(kmer_found)
			kmer_found = [0, kmer_found[1]]
		if kmer_found[0] != 0:
			if kmer_found[0] in detail_dict:
				detail_dict[kmer_found[0]] += kmer_found[1]
			else:
				detail_dict[kmer_found[0]] = kmer_found[1]
	return detail_dict


def create_list_detailed_kraken(name, status, result_kr, length, detail_dict):
	"""Columns : 
	name, status, result_kr, length_contig, nb_hits, taxon, scientific name"""
	res = []
	row_df = [name, status, result_kr, length]

	# aller chercher le nom de l'organisme avec la classe KrakenResults
	for taxon in detail_dict.keys():
		try:
			info_taxon = k_result.tax[taxon]
			#print("Found %d" %taxon)
		except KeyError:
			#print("Taxon %d not found : fetching from Ensembl" %taxon)
			info_taxon = k_result.tax.fetch_by_id(taxon)

		try:
			name_tax = info_taxon['scientific_name']
		except KeyError:
			name_tax = "Unknown"
			#print(info_taxon)
		except TypeError:
			name_tax = "Unknown"

		nb_hits = detail_dict[taxon]
		res.append(row_df + [nb_hits ,taxon, name_tax])

	return res


################################ INPUT DATA ##############################################################################################

k_result = kraken.KrakenResults(filename_kraken)

# save result without detailed info
k_result.df.to_csv(filename_output)

################################ EXECUTE ##############################################################################################

# # contient les informations par contig (dans l'ordre, pas de nom)
# k_res.df

# # df avec les numeros des taxons et le nombre de contigs trouves pour ce taxon
# k_res.taxons

# # renvoie un dictionnaire avec des infos sur un taxons
# tt_name = k_res.tax[83560]
# # renvoie le nom de l'organisme
# tt_name['scientific_name']


list_detailed_kraken = []

f = open(filename_kraken, 'r')

contig_info = f.readline()
while contig_info:
	contig_list = contig_info.replace("\n","").split("\t")
	status = contig_list[0]
	name = contig_list[1]
	result_kr = contig_list[2]
	length = contig_list[3]
	detail_kr = contig_list[4]
	# print("")
	# print(status, name, result_kr,length)

	# classified sequences : get details
	if status == "C":
		detail_kr = detail_kr.split(" ")
		detail_dict = parse_details_kraken(detail_kr)
		res_one_contig = create_list_detailed_kraken(name, status, result_kr, length, detail_dict)
	else:
		res_one_contig = [[name, status, result_kr, length, 0 , 0, "Unknown"]]
	list_detailed_kraken = list_detailed_kraken + res_one_contig
	# print(detail_dict)
	contig_info = f.readline()



f.close()


df_detailed_kraken = pd.DataFrame(list_detailed_kraken)
df_detailed_kraken.columns = ["name","status","result_kr","length_contig","nb_hits","taxon","scientific_name"]

df_detailed_kraken.to_csv(filename_output.replace(".csv","_detailed.csv"),index=None)

################################ PLOTS ##############################################################################################

#### in progress
"""

name_contig = "ENA_CP002217.fa_start_3633388_len_71166"
df_to_plot = df_detailed_kraken[df_detailed_kraken['name'] == name_contig].copy()
df_to_plot.sort_values(by="nb_hits",ascending=True,inplace=True)
x = range(df_to_plot.shape[0])
pylab.barh(x,df_to_plot["nb_hits"])
# pylab.xticks(x, df_to_plot["scientific name"], rotation='vertical')
for i, v in enumerate(df_to_plot["scientific_name"]):
    pylab.text(50, i,str(v), color='black' )
pylab.title(name_contig)
pylab.yticks([])
pylab.xlabel("Number of k-mers found")
pylab.ylabel("Species")
pylab.show()

"""


