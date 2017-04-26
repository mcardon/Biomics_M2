# -*- coding: utf-8 -*-
# arg1 : file.fastq
# arg2 : fof_databases : path/to/database (in order of use)
# arg3 : name of output directory (will be created)
# arg4 : threads


################################ IMPORT ################################################################################################
import sys
from sequana import kraken

################################ PARAMETERS ##############################################################################################

filename_fastq        = str(sys.argv[1])
fof_databases         = str(sys.argv[2])
prefix_output_results = str(sys.argv[3])
threads               = int(sys.argv[4])



kh = kraken.KrakenHierarchical(filename_fastq, fof_databases, threads=threads,output_directory=prefix_output_results)
kh.run()



