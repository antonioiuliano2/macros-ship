#script to recognize daughters of charm and other particles from the primary proton interaction (created on 13 May 2019)

import ROOT as r
import sys

#opening the file and getting the tree

def conversion2fedraunits(x,y,z):
   fedrax = x* cmtomicron + 62500
   fedray = y* cmtomicron + 49500
   fedraz = (z - 125.56649)*cmtomicron
   return fedrax, fedray, fedraz

def decaylen(mothertrack, track):
  """ Estimation of charm decay length: distance between start of mother and start of daughter """
  mx = mothertrack.GetStartX()
  my = mothertrack.GetStartY()
  mz = mothertrack.GetStartZ()

  tx = track.GetStartX()
  ty = track.GetStartY()
  tz = track.GetStartZ()

  return r.TMath.Sqrt(pow(tx-mx,2)+pow(ty-my,2)+pow(tz-mz,2))

def getdaughtertracks(inputtree,eventnumber):

 #interesting pdgcodes to check
 intermediatelist = [223, 3332, 3224, 331, 221, 20213, 3212, 213, 113, 2224, 323] #particles with very short lifetime
 signallist = [431, 411, 4122, 421, 4132, 4232, 4332, 441] #charmed hadrons

 #for recognizing particles name
 pdg = r.TDatabasePDG.Instance()
#information to be stored
 tracksID = []
 charmIDs = []
 ndaughters = {} #how many charged daughters are expected?
 decaylength = {}

 inputtree.GetEntry(eventnumber)
 mctracks = inputtree.MCTrack
#****************************** loop on tracks*******************
 for j, track in enumerate(mctracks): 
   name = "UNKNOWN"

   pdgcode = track.GetPdgCode()
   motherID = track.GetMotherId()
   charge = 0.

   if pdg.GetParticle(pdgcode):
    name = pdg.GetParticle(pdgcode).GetName()
    charge = pdg.GetParticle(pdgcode).Charge()

   if (r.TMath.Abs(pdgcode) in signallist): #charmed hadron
    charmIDs.append(j)
    ndaughters[j] = 0 #new charm, starting counter for daughters
    decaylength[j] = 0.

   cmtomicron = 1E+4

   if (motherID >= 0):

    mothertrack = mctracks[motherID]
    motherpdg = mothertrack.GetPdgCode()
    #intermediate state, continue looking for mother
    while r.TMath.Abs(motherpdg) in intermediatelist:

     motherID = mothertrack.GetMotherId()
     mothertrack = mctracks[motherID]
     motherpdg = mothertrack.GetPdgCode()

  #checking if daughter of charm
    if ((r.TMath.Abs(motherpdg) in signallist) and (r.TMath.Abs(pdgcode) not in intermediatelist)and (abs(charge)>0.)):
      #print(j, pdgcode, momentum, motherpdg, motherID, startx, starty, startz)
     tracksID.append(j)
     ndaughters[motherID] = ndaughters[motherID] + 1
     decaylength[motherID] = decaylen(mothertrack,track)
 return charmIDs, tracksID, ndaughters, decaylength

