# -*- coding: utf-8 -*-
# arg1 : time.csv
# arg2 : prefix filename for plots (if show : don't save anything, show plots)

################################ IMPORT ################################################################################################
import pylab
from sequana.lazy import pandas as pd

import sys

################################ PARAMETERS ##############################################################################################

csv_file        = str(sys.argv[1])
filename_output = str(sys.argv[2])

fontsize = 16

################################ INPUT DATA ##############################################################################################

df_time = pd.read_csv(csv_file)
print(df_time)

for align in df_time["Aligneur"].unique():
	df_to_plot = df_time[df_time["Aligneur"] == align]
	pylab.plot(df_to_plot["Threads"], df_to_plot["real_time_H"], linestyle='--', label=align)
pylab.xlabel("Threads",fontsize=fontsize)
pylab.ylabel("Temps (heures)",fontsize=fontsize)
pylab.title("Temps d'éxecution",fontsize=fontsize)
pylab.legend(fontsize=fontsize)
pylab.tight_layout()
if filename_output == "show":
	pylab.show()
else:
	pylab.savefig(filename_output, dpi=150)

pylab.clf()

