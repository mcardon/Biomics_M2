# -*- coding: utf-8 -*-
# arg1 : csv file input
# arg2 : filename output.markdown
# arg3 : optionnal : header? ( if anything here, there is a header)

################################ IMPORT ################################################################################################

import sys

################################ PARAMETERS ##############################################################################################

file_input = str(sys.argv[1])
filename_output = str(sys.argv[2])

if len(sys.argv) > 3:
	header = True

################################ FUNCTIONS ##############################################################################################


def line_to_write(list_data):
	l_out = "| " + " | ".join(list_data) + " |  \n"
	return l_out


################################ EXECUTE ##############################################################################################

f = open(file_input, 'r')
f_out = open(filename_output, 'w')

l = f.readline()
first_line = l.replace("\n","").split(",")

if header:
	# first line contains header : write header
	f_out.write(line_to_write(first_line))
	# go to next line (data)
	l = f.readline()

else:
	# first line contains data : write header first
	num_head = [str(i) for i in range(len(first_line))]
	f_out.write(line_to_write(num_head))

# write separator
f_out.write(line_to_write(["---"]*len(first_line)))

while(l):
	f_out.write(line_to_write(l.replace("\n","").split(",")))
	l = f.readline()

f.close()
f_out.close()



