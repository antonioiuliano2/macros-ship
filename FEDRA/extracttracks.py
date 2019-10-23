import fedrarootlogon
import sys
import ROOT

inputfile = ROOT.TFile.Open(sys.argv[1],"UPDATE")
inputtree = inputfile.Get("tracks")

inputtree.SetAlias("t","t.")

nentries = inputtree.GetEntries()
tracks = ROOT.TClonesArray("EdbSegP") #try with a TClonesArray
tracks_COV = ROOT.TClonesArray("TMatrixD")

copiedseg = ROOT.EdbSegP()
inputtree.Branch("track",tracks) 

inputtree.Branch("track_COV",tracks_COV)
#we add a branch with the copied track as EdbSegP, to solve reading issue
print ("Starting loop over {} events".format(nentries))
for ientry in range(nentries):
 #clearing TClonesArrays
 tracks.Clear("C")
 tracks_COV.Clear("C")


 copiedseg = ROOT.EdbSegP()
 inputtree.GetEntry(ientry)
 track = inputtree.t
 #getting covariance matrix and adding it to the TClonesArray
 covtrack = track.COV()
 tracks_COV[0] = covtrack
 #copying track branch into a TClonesArray of EdbSegP
 #copiedseg.Copy(track)
 copiedseg = ROOT.EdbSegP(track)
 tracks[0] = copiedseg

 inputtree.Fill()
 if ((ientry % 100000) == 0) :
  print ("arrived at entry ",ientry)
print ("Finished loop over {} events".format(nentries))
inputfile.Write()
inputfile.Close()
