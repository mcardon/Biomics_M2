# -*- coding: utf-8 -*-
# arg1 : bam file
# arg2 : filename output.bam
# arg3 : stride

################################ IMPORT ################################################################################################
import os
import sys
#from sequana import *

scriptpath = "/home/mcardon/Mel/Code/sequana/sequana/pacbio.py"
sys.path.append(os.path.abspath(scriptpath))
import pacbio

################################ PARAMETERS ##############################################################################################

filename_BAM = str(sys.argv[1])
filename_output = str(sys.argv[2])
stride_nb = int(sys.argv[3])


################################ INPUT DATA ##############################################################################################

raw_BAM = pacbio.PacBioInputBAM(filename_BAM)
raw_BAM.stride(filename_output, stride=stride_nb)

################################ EXECUTE ##############################################################################################


################################ PLOTS ##############################################################################################









