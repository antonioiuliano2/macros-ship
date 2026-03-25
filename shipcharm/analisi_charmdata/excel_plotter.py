# script to plot information stored in a excel data sheet, after importing them to panda
from __future__ import division # for integer division -> python 2 will lead to an integer

import os
import pandas as pd
from matplotlib import pyplot as plt

homepath = os.environ['HOME']
#converting it into a panda dataframe
df = pd.read_excel(homepath+"/Dropbox/Charmdata/alignment_worklist.xlsx",'Summary_forexport')

fig, ax = plt.subplots()

df.index = df['RUN'] #index is displayed as plot x axis, automatically!

percdone = df['Ndone']/df['Npairs'] * 100
percgood = df['Ngood']/df['Npairs'] * 100

plt.xlabel('RUN')
plt.ylabel('Percentage (%)')

ax.plot(percdone,'b.',label='scanned and processed for alignment')
ax.plot(percgood,'g.',label='aligned successfully')

plt.legend(fontsize='x-large')
plt.show()
