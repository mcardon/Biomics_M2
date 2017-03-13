# -*- coding: utf-8 -*-
# arg1 : Read Length file (RL part of output from samtools stats)

#################################
#         IMPORT
#################################

import pandas as pd
import matplotlib.pyplot as plt
import sys

filename = str(sys.argv[1])
print(filename)

data = pd.read_csv(filename, sep="\t", header=None)
data.columns = ["dummy","length", "count"]

plt.plot(data["length"], data["count"])
plt.ylim(0,100)
plt.xlim(0,50000)
plt.xlabel("Read length")
plt.ylabel("Nb of reads")
plt.title("Number of reads per length")
plt.savefig("Nb_reads_len_scraps_cut100.png")
plt.show()


