#script to cut cbmsim tree according to given selection
import ROOT as r
import numpy as np

nfiles = 1000
#selection
xmin = -20.00 #-20.05
xmax = 20.00 #20.05

ymin = -20.00 #-20.05
ymax = 20.00 #20.05

simchain = r.TChain("cbmsim")

for ifile in range(nfiles):
 simchain.Add("{}/ship.conical.Genie-TGeant4.root".format(ifile+1)) #file names start from 1
#cloning the tree, without any entry
copyfile = r.TFile.Open(("inECC_ship.conical.Genie-TGeant4.root"),"RECREATE")
copytree = simchain.CloneTree(0)
#adding information (I WANT TO BE ABLE TO GET THE ORIGINAL GENIE EVENT)

GenieEventID = np.zeros(1,dtype=np.intc) # number of tree (aka condor job + 1)
copytree.Branch("GenieEventID",GenieEventID,"GenieEventID/I")


nevents = simchain.GetEntries()
print("Start processing ",nevents,"in region x in ",xmin,xmax," and y in ",ymin, ymax)
#loop into events
for ievent in range(nevents):
 simchain.GetEvent(ievent)
 tracks = simchain.MCTrack
 startx = tracks[0].GetStartX()
 starty = tracks[0].GetStartY()
 #adding my additional variables
 GenieEventID[0] = simchain.GetTreeNumber() + 1

 if (startx >= xmin) and (startx <= xmax):
  if (starty >= ymin) and (starty <= ymax):
   copytree.Fill()

copytree.AutoSave()
copyfile.Close()
