# -*- coding: utf-8 -*-
# arg1 : file.fastq
# arg2 : fof_databases : path/to/database (in order of use)
# arg3 : name of output directory (will be created)
# arg4 : threads


################################ IMPORT ################################################################################################
import os
import sys
from sequana import kraken
from easydev import DevTools, execute, TempFile, md5
################################ PARAMETERS ##############################################################################################

filename_fastq        = str(sys.argv[1])
fof_databases         = str(sys.argv[2])
prefix_output_results = str(sys.argv[3])
threads               = int(sys.argv[4])

################################ EXECUTE ##############################################################################################


class KrakenHierarchical(object):
    """
    This runs Kraken on a fastq file with multiple k-mer databases in a hierarchical way.
    Unclassified sequences with the first database are input for the second, and so on.


    """
    def __init__(self, filename_fastq, fof_databases, threads=1,output_directory="./kraken_hierarchical/", keep_temp_files=False):
        """.. rubric:: **constructor**

        :param filename_fastq: fastq file to analyse
        :param fof_databases: file with absolute paths to databases (one per line, in order of desired use)
        :param threads: number of threads to be used by Kraken
        :param output_directory: name of the output directory
        :param keep_temp_files: bool, if True, will keep intermediate files from each Kraken analysis, and save html report at each step

        """
        self.filename_fastq = filename_fastq
        with open(fof_database, 'r') as fof:
            self.databases = [absolute_path.split('\n')[0] for absolute_path in fof.readlines()]
        self.threads = threads
        self.output_directory = output_directory
        self.keep_temp_files = keep_temp_files

        os.mkdir(output_directory)

        # list of input fastq files
        self._list_file_fastq = [filename_fastq]
        # list of all output to merge at the end
        self._list_kraken_output = []

    def _run_one_analysis(self, i, last_analysis=False):
        """
        """
        analysis = KrakenAnalysis(self._list_file_fastq[i], self.databases[i], self.threads)

        file_kraken_class  = self.output_directory + "kraken_%d.out" %(i)

        if last_analysis:
            analysis.run(output_filename=file_kraken_class)
        else:
            file_fastq_unclass = self.output_directory + "unclassified_%d.fastq" %(i)
            analysis.run(output_filename=file_kraken_class, output_filename_unclassified=file_fastq_unclass)
            self._list_file_fastq.append(file_fastq_unclass)
            self._list_kraken_output.append(file_kraken_class)

        if self.keep_temp_files:
            result = kraken.KrakenResults(file_kraken_class)
            result.to_js("%skrona_%d.html" %(prefix_output_results,i))

    def run(self):
        """
        """
        # Iteration over the databases
        for i in range(len(self.databases)):
            is_last = (i == (len(self.databases)-1) )
            self._run_one_analysis(i, is_last)

        # concatenate all files
        file_output_final = self.output_directory + "kraken_final.out"
        with open(file_output_final, 'w') as outfile:
            for fname in self._list_kraken_output:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        # create html report
        result = kraken.KrakenResults(file_output_final)
        result.to_js("%skrona_final.html" %self.output_directory)

        if not self.keep_temp_files:
            # remove kraken intermediate files
            for f_temp in self._list_kraken_output:
                os.remove(f_temp)
            # remove unclassified files
            for f_temp in self._list_file_fastq[1:len(self._list_file_fastq)]:
                os.remove(f_temp)


