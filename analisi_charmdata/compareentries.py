import ROOT as r
import sys,os
import matplotlib.pyplot as plt
import pandas as pd

runnumber = sys.argv[1]

spillcodes="/afs/cern.ch/work/a/aiuliano/public/Charmdata/spillcodes/spillcodes_run"+runnumber+".csv"
datadir="/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rootdata/RUN_8000_"+runnumber+"/"
scifidatadir="/eos/experiment/ship/data/charmxsec/DATA_0900/rootdata/RUN_0900_"+runnumber+"/"

#nentriesscifi = []

counters = []

nentries1 = []
nentries2 = []
ndeltaentries = []
spillcodelist= []

with open(spillcodes) as f:
 for i, spillcode in enumerate(f):
  counters.append(i)

  spillcode1 = spillcode.strip() #removes \n, in order to make it usable
  spillcodelist.append(str(spillcode1)) 

  file1loc = "ls "+ datadir+"*"+str(spillcode1)+"*.root"  
  FILE1 = os.popen(file1loc).read()
  FILE1 = FILE1.strip()

  file2loc = "ls "+ scifidatadir+"*"+spillcode1+"*.root"  
  FILE2 = os.popen(file2loc).read()
  FILE2 = FILE2.strip()
  #getting number of entries from trees
  file1 = r.TFile.Open(FILE1)
  tree1 = file1.Get("cbmsim")   

  file2 = r.TFile.Open(FILE2)
  tree2 = file2.Get("cbmsim")
 # nentriesscifi.append(tree2.GetEntries())
  nentries1.append(tree1.GetEntries())
  nentries2.append(tree2.GetEntries())
  ndeltaentries.append(tree1.GetEntries() - tree2.GetEntries())

  file1.Close();
  file2.Close();

table = {'spillcodelist':spillcodelist,'nentries1':nentries1,'nentries2':nentries2,'ndeltaentries':ndeltaentries}

entriestable = pd.DataFrame(table,columns = ['spillcodelist','nentries1','nentries2','ndeltaentries'])
entriestable.to_csv('entriestable'+runnumber+'.csv')

#plt.plot(nentriesscifi, "r*", label = "DATA_900(SciFi)")
plt.plot(ndeltaentries, "b*", label = "DATA_8000(Others)- DATA_900(SciFi)")
plt.xticks(counters, spillcodelist, rotation = "vertical")
plt.ylabel("nentries_others - nentries_scifi")
plt.yscale("log") #for log scale
plt.title("Comparison of entries for run "+runnumber)
#plt.legend()
plt.show()

