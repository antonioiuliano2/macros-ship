#script to recognize daughters of charm and other particles from the primary proton interaction (created on 13 May 2019)

import ROOT as r
from rootUtils import bookHist
import sys

#opening the file and getting the tree

def getdaughtertracks(inputtree,eventnumber):

 #interesting pdgcodes to check
 intermediatelist = [223, 3332, 3224, 331, 221, 20213, 3212, 213, 113, 2224, 323] #particles with very short lifetime
 signallist = [431, 411, 4122, 421, 4132, 4232, 4332, 441] #charmed hadrons

 #for recognizing particles name
 pdg = r.TDatabasePDG.Instance()

 tracksID = []

 inputtree.GetEntry(eventnumber)
 mctracks = inputtree.MCTrack
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

   cmtomicron = 1E+4
   startx = track.GetStartX()* cmtomicron + 62500;
   starty = track.GetStartY()* cmtomicron + 49500;
   startz = (track.GetStartZ() - 125.56649)*cmtomicron;

   if (motherID >= 0):

    mothertrack = mctracks[motherID]
    motherpdg = mothertrack.GetPdgCode()
    #intermediate state, continue looking for mother
    while r.TMath.Abs(motherpdg) in intermediatelist:

     motherID = mothertrack.GetMotherId()
     mothertrack = mctracks[motherID]
     motherpdg = mothertrack.GetPdgCode()

  #checking if daughter of charm
    if ((r.TMath.Abs(motherpdg) in signallist) and (r.TMath.Abs(pdgcode) not in intermediatelist)):
     print j, pdgcode, momentum, motherpdg, motherID, startx, starty, startz
     tracksID.append(j)
 return tracksID

