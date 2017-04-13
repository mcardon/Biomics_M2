# -*- coding: utf-8 -*-
# arg1 : fof of all genome rotations (no path, only filenames)
# arg2 : filename fasta file to evaluate (contigs.fasta)
# arg3 : output prefix for script (no extension)

### rotations of the genome must be in ./../ref_rotate/ directory (change this if necessary)

################################ IMPORT ################################################################################################

import time
import sys

################################ PARAMETERS ############################################################################################

fof_ref = str(sys.argv[1])
filename_input = str(sys.argv[2])
filename_output = str(sys.argv[3])

directory_ref_rotate = "./../ref_rotate/"

################################ IMPORT ################################################################################################


# list of all files
f = open(fof_ref, 'r')
list_files = []
for name in f.readlines():
	f_name = name.split('\n')[0]
	list_files.append(f_name)
f.close()


today = time.strftime("%Y_%m_%d")
header_script = ["#!/bin/sh","#SBATCH --qos=fast","#SBATCH -c 2","#SBATCH --mem-per-cpu=10000","#SBATCH -p dedicated"]
header_script.append("#SBATCH -o ./logs/%s_%s_ref_rotate.out -e ./logs/%s_%s_ref_rotate.err" % (today, filename_output.split("/")[-1], today, filename_output.split("/")[-1]))

body_script = []
for i in range(len(list_files)):
	body_script.append("srun blasr --nproc 2 --minMatch 15 --maxMatch 20 --advanceHalf --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut %s %s%s -m 5 --out %s%s.blasr5.out || exit %d" % (filename_input,directory_ref_rotate,list_files[i], filename_input.replace(".fasta",""), list_files[i].replace(".fasta",""), i+1))


tail_script = ["exit 0"]


script = header_script + body_script + tail_script

f_out = "%s_%s.sh" %(filename_output, today)

f = open(f_out, 'w')
f.write("\n".join(script))
f.close()

