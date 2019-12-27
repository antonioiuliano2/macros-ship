#script to cut cbmsim tree according to given selection
import ROOT as r
import sys

#selection
zmin = 121.8365
zmax = 125.5500

simfile = r.TFile.Open(sys.argv[1])
simtree = simfile.Get("cbmsim")
#cloning the tree, without any entry
copyfile = r.TFile.Open(("inECC_"+sys.argv[1]),"RECREATE")
copytree = simtree.CloneTree(0)

nevents = simtree.GetEntries()
#loop into events
for ievent in range(nevents):
 simtree.GetEvent(ievent)
 tracks = simtree.MCTrack
 startz = tracks[0].GetStartZ()

 if (startz >= zmin) and (startz <= zmax):
  copytree.Fill()

copytree.AutoSave()
copyfile.Close()
