# -*- coding: utf-8 -*-
# arg1 : file.fastq
# arg2 : fof_databases : path/to/database (in order of use)
# arg3 : name of output directory (will be created)
# arg4 : threads


################################ IMPORT ################################################################################################
import os
import sys
from sequana import kraken

################################ PARAMETERS ##############################################################################################

file_fastq            = str(sys.argv[1])
fof_database          = str(sys.argv[2])
prefix_output_results = str(sys.argv[3])
threads               = int(sys.argv[4])

################################ EXECUTE ##############################################################################################

# fastq files
list_file_fastq = [file_fastq]
# list of all databases paths
f = open(fof_database, 'r')
list_databases = [name.split('\n')[0] for name in f.readlines()]
f.close()
# list of all output to merge at the end
list_kraken_output = []

# create output directory
os.mkdir(prefix_output_results)

# Iteration over the N-1 databases
for i in range(len(list_databases)-1):
	analysis = kraken.KrakenAnalysis(list_file_fastq[i], list_databases[i], threads)

	file_kraken_class  = prefix_output_results + "kraken_%d.out" %(i)
	file_fastq_unclass = prefix_output_results + "unclassified_%d.fastq" %(i)
	analysis.run(output_filename=file_kraken_class, output_filename_unclassified=file_fastq_unclass)
	list_file_fastq.append(file_fastq_unclass)
	# XXX remove for final : save krona.html
	result = kraken.KrakenResults(file_kraken_class)
	result.to_js("%skrona_%d.html" %(prefix_output_results,i))



# Last database
i=len(list_databases)-1
# Run with last DB and output everything
analysis = kraken.KrakenAnalysis(list_file_fastq[i], list_databases[i], threads)
file_kraken_class = prefix_output_results + "kraken_%d.out" %(i)
analysis.run(output_filename=file_kraken_class)
# XXX remove for final : save krona.html
result = kraken.KrakenResults(file_kraken_class)
result.to_js("%skrona_%d.html" %(prefix_output_results,i))

# concatenate all files
file_output_final = prefix_output_results + "kraken_final.out"
filenames_kraken  = [prefix_output_results + "kraken_%d.out" %(i) for i in range(len(list_databases))]
with open(file_output_final, 'w') as outfile:
    for fname in filenames_kraken:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

# create html report
result = kraken.KrakenResults(file_output_final)
result.to_js("%skrona_final.html" %prefix_output_results)

# for final version of code, remove temporary kraken files and unclassified.fastq



