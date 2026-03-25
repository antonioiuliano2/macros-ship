#Testing branch uprooting, to load them into keras, created on 26 September 2019
import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

signaltree = uproot.open("vertices_MC_modified.root")["vtx"]
backgroundtree = uproot.open("vertextree_mateisim.root")["vtx"]

def builddataframe(mytree):
 '''building the Pandas dataframe from the ROOT tree'''
 probability = mytree.array("probability")
 molteplicity = mytree.array("n")

 vx = mytree.array("vx")
 vy = mytree.array("vy")
 vz = mytree.array("vz")

 if "vtx_max_aperture" in mytree.keys():
  aperture = mytree.array("vtx_max_aperture")
 else: 
  aperture = mytree.array("maxaperture")

 data = pd.DataFrame({"Probability":probability,"Vx":vx,"Vy":vy,"Vz":vz,"MaxAperture":aperture},columns = ['Probability','Vx','Vy','Vz','MaxAperture'])
 
 return data

signaldata = builddataframe(signaltree)
backgrounddata = builddataframe(backgroundtree)

fig,ax = plt.subplots()
hnp = ax.hist(signaldata['Probability'],10,[0,1])
ax.set_xlabel("probability")

#applying a filter
probabilitycut = 0.9
signaldata = signaldata[signaldata['Probability']>probabilitycut]
backgrounddata = backgrounddata[backgrounddata['Probability']>probabilitycut]

#histogram drawing
fig1,ax1 = plt.subplots()
htxty = ax1.hist2d(signaldata["Vx"],signaldata["Vy"],[120,100],[[0,120E+3],[0,100E+3]])
fig1.colorbar(htxty[3])
ax1.set_xlabel("vx[micron]")
ax1.set_ylabel("vy[micron]")

fig2,ax2 = plt.subplots()
hvz = ax2.hist(signaldata["Vz"])
ax2.set_xlabel("vz[micron]")


#data normalization: mean 0 and var 1


def norm(x):
  x.pop('Probability') #already used for the first cut
  stats = x.describe().transpose()
  return ((x - stats['mean'])/stats['std'])
signalstats=signaldata.describe().transpose()
backgroundstats=backgrounddata.describe().transpose()

normed_signaldata = norm(signaldata) #check with describe
normed_backgrounddata = norm(backgrounddata)

normed_signalstats=normed_signaldata.describe().transpose()
normed_backgroundstats=normed_backgrounddata.describe().transpose()
#plotting output to excel for report
excelwriter = pd.ExcelWriter("root2dataframestats.xls",engine = 'xlsxwriter')
signalstats.to_excel(excelwriter,sheet_name='Signal')
backgroundstats.to_excel(excelwriter,sheet_name='Background')
normed_signalstats.to_excel(excelwriter,sheet_name='NormedSignal')
normed_backgroundstats.to_excel(excelwriter,sheet_name='NormedBackground')

#comparisons between maxaperture and maxip
fig3,ax3 = plt.subplots()
ax3.hist(normed_signaldata['MaxAperture'],alpha=0.5,label='signalsim')
ax3.hist(normed_backgrounddata['MaxAperture'],alpha=0.5,label='backgroundsim')
ax3.set_xlabel('Norm_MaxAperture')

ax3.legend()

excelwriter.save()
plt.show()


