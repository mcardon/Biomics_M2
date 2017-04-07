# -*- coding: utf-8 -*-
# arg1 : bam file
# arg2 : filename output.bam
# arg3 : threshold

################################ IMPORT ################################################################################################
import os
import sys
#from sequana import *

scriptpath = "/home/mcardon/Mel/Code/sequana/sequana/"
sys.path.append(os.path.abspath(scriptpath))
import pacbio

################################ PARAMETERS ##############################################################################################

filename_BAM = str(sys.argv[1])
filename_output = str(sys.argv[2])
threshold = int(sys.argv[3])

################################ INPUT DATA ##############################################################################################

raw_BAM = pacbio.PacBioInputBAM(filename_BAM)


################################ EXECUTE ##############################################################################################

raw_BAM.filter_length(filename_output, threshold)


################################ PLOTS ##############################################################################################









