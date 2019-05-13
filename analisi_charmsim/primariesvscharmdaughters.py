#script to recognize daughters of charm and other particles from the primary proton interaction (created on 13 May 2019)

import ROOT as r
from rootUtils import bookHist
import sys

#opening the file and getting the tree

filename = sys.argv[1]

inputfile = r.TFile.Open(filename)
tree = inputfile.Get("cbmsim")

nevents = tree.GetEntries()
#interesting pdgcodes to check
intermediatelist = [223, 3332, 3224, 331, 221, 20213, 3212, 213, 113, 2224, 323] #particles with very short lifetime
signallist = [431, 411, 4122, 421, 4132, 4232, 4332, 441] #charmed hadrons

#label for the outputfile
print "IEVENT, ITRACK, PDGCODE,MOMENTUM, MOTHERPDG, MOTHERID"

#for recognizing particles name
pdg = r.TDatabasePDG.Instance()

histos = {}
bookHist(histos, "hlen", "Decay length", 30, 0, 3)

#******************************loop on events*********************
for i in range(100): 
 tree.GetEntry(i)
 mctracks = tree.MCTrack
#****************************** loop on tracks*******************
 for j, track in enumerate(mctracks): 
  name = "UNKNOWN"

  pdgcode = track.GetPdgCode()
  px = track.GetPx()
  py = track.GetPy()
  pz = track.GetPz()
  momentum = track.GetP()
  motherID = track.GetMotherId()

  if pdg.GetParticle(pdgcode):
   name = pdg.GetParticle(pdgcode).GetName()

  startz = track.GetStartZ()

  if (motherID == -1): #daughter of primary proton
   print i, j, pdgcode, momentum, 2212, -1 

  elif (motherID > 0):

   mothertrack = mctracks[motherID]
   motherpdg = mothertrack.GetPdgCode()
   #intermediate state, continue looking for mother
   while r.TMath.Abs(motherpdg) in intermediatelist:

     motherID = mothertrack.GetMotherId()
     mothertrack = mctracks[motherID]
     motherpdg = mothertrack.GetPdgCode()

  #checking if daughter of charm
   if ((motherpdg in signallist) and (pdgcode not in intermediatelist)):
    print i, j, pdgcode, momentum, motherpdg, motherID

    decaylen = pow(pow(track.GetStartX() - mothertrack.GetStartX(),2)+pow(track.GetStartY() - mothertrack.GetStartY(),2)+pow(track.GetStartZ() - mothertrack.GetStartZ(),2),0.5)
    histos["hlen"].Fill(decaylen)

histos["hlen"].Draw()
