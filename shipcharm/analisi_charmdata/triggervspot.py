import matplotlib.pyplot as plt
import numpy as np
import ROOT
import os #we can get list of file with listdir!
'''Script to compare entries from SHiP-charm data and expected pots'''
inputfile = ROOT.TFile.Open("/eos/experiment/ship/data/charmxsec/bookkeeping/charm_spills_predictions.root") #file to take numberss of pots and spills from

inputtree = inputfile.Get("spill")
#writing a log
outputfile = open("/eos/experiment/ship/data/charmxsec/bookkeeping/triggervspot.csv",'w')

outputfile.write("nrun,runname,spillcode,ntriggers,nPOTs\n")
#loop on runs
for run in inputtree:
 nrun = run.runcode
 runname = run.name
 datapath = "/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rootdata/RUN_8000_{}/".format(nrun)
#liste di POT e di spillcode estratte dal tree

#getting pots and spillcodes

 pots = run.pot
 spillcodes = run.spillcode

 POTslist = []
 spillcodelist = []
 filenames = []
 filelist = os.listdir(datapath) #list of spill for a given run

 entrieslist = []
 for pot,spillcode in zip(pots,spillcodes):

  if (pot > 0):
   POTslist.append(pot)
   for name in filelist: #finding the spill file corresponding to spillcode
    if name.__contains__(str(spillcode)):
     filenames.append(name)
     spillcodelist.append(spillcode)
 print ("Getting entries for run {}".format(nrun)) 
 for ispill, filename in enumerate(filenames):

  datafile = ROOT.TFile.Open(datapath+filename)
  datatree = datafile.Get("cbmsim")
  
  entries = datatree.GetEntries()
  entrieslist.append(entries)
  #are the triggers more than the pots?
  outputfile.write("{},{},{},{},{}\n".format(nrun,runname,spillcodes[ispill],entries,POTslist[ispill]))
  #if entries > POTslist[ispill]:
  #	print ("Run {}, spill {}: we have {} triggers vs {} POTs".
  #         format(nrun,spillcodes[ispill], entries, POTslist[ispill]))
 #plot to show comparison for each spill
 figure0 = plt.figure()
 plt.title("Comparing entries and POTs for run {}, Charm {}".format(nrun,runname))
 plt.plot(entrieslist,"b*",label="Nentries_{}".format(nrun))
 plt.plot(POTslist,"ro",label="nPOTs_{}".format(nrun))
 plt.xticks(np.arange(len(POTslist)),  spillcodelist, rotation='vertical',fontsize="x-small")
 plt.ylabel("POTs")
 plt.legend()
 plt.savefig('/eos/user/a/aiuliano/Synched/Archivio_cronologico/Novembre 2019/entriesandpots_{}.png'.format(nrun))
 plt.close()
#end of run loop, close output file
outputfile.close()
