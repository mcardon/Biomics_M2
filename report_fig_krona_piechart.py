# -*- coding: utf-8 -*-
# arg1 : data.csv
# arg2 : prefix filename for plots (if show : don't save anything, show plots)

################################ IMPORT ################################################################################################
import pylab
from sequana.lazy import pandas as pd

import sys

################################ PARAMETERS ##############################################################################################

csv_file        = str(sys.argv[1])
filename_output = str(sys.argv[2])

fontsize = 16
pylab.rcParams['font.size'] = 8
nb_row = 2
nb_col = 5

colors = ["#6fcfcf","#cf9f6f","#cf6f6f"]

################################ PARAMETERS ##############################################################################################

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%\n({v:d})'.format(p=pct,v=val)
    return my_autopct



################################ INPUT DATA ##############################################################################################

df_data = pd.read_csv(csv_file)
# print(df_data)

fig1, axarr = pylab.subplots(nb_row,nb_col, figsize=(nb_col*2, nb_row*2))

i=0
for (db,tp) in df_data.loc[:,("database","type")].drop_duplicates().values:
	print(db,tp)
	df_to_plot = df_data[(df_data["database"] == db) & (df_data["type"] == tp) ]
	ax = axarr[int(i/nb_col)][i%nb_col]
	sizes = df_to_plot["nb_reads"]
	labels = df_to_plot["result"]
	explode = (0.1,0.,0.)
	total = sizes.sum()
	ax.pie(sizes, explode=explode, colors=colors, autopct=make_autopct(sizes),shadow=False, startangle=-45)
	ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	i += 1
	# pylab.plot(df_to_plot["Threads"], df_to_plot["real_time_H"], linestyle='--', label=align)
if filename_output == "show":
	fig1.show()
else:
	fig1.savefig(filename_output,dpi=300)
	fig1.clf()






