# -*- coding: utf-8 -*-
# arg1 : low_coverage.csv (output of sequana_coverage standalone)
# arg2 : filename for output plot (if "show" dont save, show instead)

################################ IMPORT ################################################################################################
import pylab
from sequana.lazy import pandas as pd
import collections
import sys

################################ PARAMETERS ##############################################################################################

f_input         = str(sys.argv[1])
filename_output = str(sys.argv[2])

save_output = True
if filename_output == "show":
	save_output = False

colors = ["y","c","b","g","r","m","c"]
# colors = ["k","r","g","y","m","c"]

alpha = 0.6
nb_bars = 15
figsize = (14,8)
max_chars = 60
fontsize = 24

if save_output:
	title = filename_output.replace(".png","").replace("_"," ").replace(".txt","")
else:
	title = f_input.replace("fof_","").replace(".txt","")

################################ INPUT DATA ##############################################################################################

df = pd.read_csv(f_input)

################################ EXECUTE ##############################################################################################

names_prod = collections.Counter(df.loc[:,"product"])

list_names = []
list_count = []
for k in names_prod.keys():
	if (str(k) != "nan") & (k != "None") :
		list_names.append(str(k))
		list_count.append(int(names_prod[k]))
	# print(names_prod[k])


list_names_sort = [n for (c,n) in sorted(zip(list_count,list_names),reverse=True)][0:nb_bars]
list_count_sort = [c for (c,n) in sorted(zip(list_count,list_names),reverse=True)][0:nb_bars]


################################ PLOTS ##############################################################################################

x = range(nb_bars)

fig, ax = pylab.subplots(1,1, figsize=figsize)
pylab.barh(x,list_count_sort[::-1],color="c",alpha=alpha)
for i, v in enumerate(list_names_sort[::-1]):
	if len(v) > max_chars:
		v = v[0:max_chars] + "[...]"
	pylab.text(5, i,v, color='black',fontsize=fontsize)
pylab.title(title,fontsize=fontsize)
pylab.xlabel("Number of regions",fontsize=fontsize)
pylab.ylabel("Genbank Annotation",fontsize=fontsize)
pylab.yticks([])
pylab.tight_layout()
if save_output:
	pylab.savefig(filename_output,dpi=182)
	pylab.clf()
else:
	pylab.show()




