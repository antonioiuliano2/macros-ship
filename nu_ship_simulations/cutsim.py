#script to cut cbmsim tree according to given selection
import ROOT as r
import numpy as np

#selection
xmin = -20.00 #-20.05
xmax = 20.00 #20.05

ymin = -20.00 #-20.05
ymax = 20.00 #20.05

simfile = r.TFile.Open("sim_all.root")
simtree = simfile.Get("cbmsim")
 
#cloning the tree, without any entry
copyfile = r.TFile.Open(("inECC_ship.conical.Genie-TGeant4.root"),"RECREATE")
copytree = simtree.CloneTree(0)
#adding information (I WANT TO BE ABLE TO GET THE ORIGINAL GENIE EVENT)

GenieEventID = np.zeros(1,dtype=np.intc) # number of tree (aka condor job + 1)
copytree.Branch("GenieEventID",GenieEventID,"GenieEventID/I")


nevents = simtree.GetEntries()
print("Start processing ",nevents,"in region x in ",xmin,xmax," and y in ",ymin, ymax)
#loop into events
for ievent in range(nevents):
 simtree.GetEvent(ievent)
 tracks = simtree.MCTrack
 startx = tracks[0].GetStartX()
 starty = tracks[0].GetStartY()
 #adding my additional variables
 GenieEventID[0] = simtree.GetTreeNumber() + 1

 if (startx >= xmin) and (startx <= xmax):
  if (starty >= ymin) and (starty <= ymax):
   copytree.Fill()

copytree.AutoSave()
copyfile.Close()
