# -*- coding: utf-8 -*-
# arg1 : file.txt
# arg2 : number of output files

# splits a file in N other files, line by line

################################ IMPORT ################################################################################################
import sys

################################ PARAMETERS ##############################################################################################

filename_input = str(sys.argv[1])
N              = int(sys.argv[2])


################################ EXECUTE ##############################################################################################

filenames_output = [ filename_input.replace(".txt","_%d_%d.txt" %(i+1,N)) for i in range(N)]
f_outputs = []
# open output files
for f_out in filenames_output:
	f_outputs.append(open(f_out, 'w'))

# open input
f = open(filename_input, 'r')
l = f.readline()
i = 0

# write outputs
while l:
	f_outputs[i].write(l)
	i = (i+1)%N
	l = f.readline()

# close input
f.close()


# close output files
for f_out in f_outputs:
	f_out.close()