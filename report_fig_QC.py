# -*- coding: utf-8 -*-
# arg1 : filename genome (fasta format)
# arg2 : prefix filename for plots (if show : don't save anything, show plots)

################################ IMPORT ################################################################################################
import pylab
import numpy as np
from sequana.sequence import DNA, Repeats

import sys

################################ PARAMETERS ##############################################################################################

filename_fasta = str(sys.argv[1])
filename_output = str(sys.argv[2])

save_output = True
if filename_output == "show":
	save_output = False
