#script to recognize daughters of charm and other particles from the primary proton interaction (created on 13 May 2019)

import ROOT as r
from rootUtils import bookHist
import sys


#opening the file and getting the tree

hipcharm = r.TH1D("hipcharm","Impact parameter charm daughters",100,0,1000)
hdecaylen= r.TH1D("hdecaylen","Decay length charmed hadrons",100,0,10)
hkink= r.TH1D("hkink","Kink angle charm and daughters",100,0,1)

def decaylen(parenttrack, daughtertrack):

 daughterstartx = daughtertrack.GetStartX()
 daughterstarty = daughtertrack.GetStartZ()
 daughterstartz = daughtertrack.GetStartY()


 parentstartx = parenttrack.GetStartX()
 parentstarty = parenttrack.GetStartZ()
 parentstartz = parenttrack.GetStartY()

 flightlength= pow(daughterstartx-parentstartx,2)+ pow(daughterstarty-parentstarty,2)+ pow(daughterstartz-parentstartz,2)
 
 return flightlength
def Kinkangle(parenttrack, daughtertrackslist):

 kink=0.

 parenttx = parenttrack.GetPx()/parenttrack.GetPz()
 parentty = parenttrack.GetPy()/parenttrack.GetPz()
 #computing the average kink angle between mother and track particles
 for daughtertrack in daughtertrackslist:

  daughtertx = daughtertrack.GetPx()/daughtertrack.GetPz()
  daughterty = daughtertrack.GetPy()/daughtertrack.GetPz()
  kink += r.TMath.Sqrt(pow(parenttx - daughtertx,2) + pow(parentty - daughterty,2))
 
 ndaughters = len(daughtertrackslist)
 return kink/ndaughters

def IPtoVertex(vertexpos, track):
 #getting vertex and start of track position as two TVector3 objects
 
 trackstartpos = r.TVector3(track.GetStartX(), track.GetStartY(), track.GetStartZ())

 px = track.GetPx()
 py = track.GetPy()
 pz = track.GetPz()
 momentum = track.GetP()

 trackdir = r.TVector3(px,py,pz)
 #now we can start Elena's loops
 delta = 0.
 for i in range(3):

  delta += (vertexpos(i) - trackstartpos(i))*trackdir(i)/momentum

 ip = 0.
 
 for i in range(3):
 
  ip = pow((vertexpos(i) - trackstartpos(i) - delta * trackdir(i)/momentum),2)
 

 return r.TMath.Sqrt(ip)


def getdaughtertracks(inputtree,eventnumber):

 #interesting pdgcodes to check
 intermediatelist = [223, 3332, 3224, 331, 221, 20213, 3212, 213, 113, 2224, 323] #particles with very short lifetime
 signallist = [431, 411, 4122, 421, 4132, 4232, 4332, 441] #charmed hadrons

 #for recognizing particles name
 pdg = r.TDatabasePDG.Instance()

 tracksID = []
 charmhadrons = []
 charmpdgs = []
 charmdaughters = [[],[]]

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
   if (r.TMath.Abs(pdgcode) in signallist): #charmed hadron
    charmhadrons.append(track)
    charmpdgs.append(pdgcode)

   if pdg.GetParticle(pdgcode):
    name = pdg.GetParticle(pdgcode).GetName()

   cmtomicron = 1E+4
   startx = track.GetStartX()
   starty = track.GetStartY()
   startz = track.GetStartZ()
   startpos = r.TVector3(startx,starty,startz)

   if (motherID == -1):
    vertex = startpos #saving position of primary vertex

   if (motherID >= 0):

    mothertrack = mctracks[motherID]
    motherpdg = mothertrack.GetPdgCode()
    #intermediate state, continue looking for mother
    while r.TMath.Abs(motherpdg) in intermediatelist:

     motherID = mothertrack.GetMotherId()
     mothertrack = mctracks[motherID]
     motherpdg = mothertrack.GetPdgCode()

  #checking if daughter of charm
    if ((motherpdg in charmpdgs) and (r.TMath.Abs(pdgcode) not in intermediatelist)):
     if momentum > 0.1: #energy cut
  
      charmindex = charmpdgs.index(motherpdg) #what charm was found?
      charmdaughters[charmindex].append(track)
      impactparameter = IPtoVertex(vertex,track)
      hipcharm.Fill(impactparameter*cmtomicron)

      tracksID.append(j)

 #loop on the two groups of charmdaughters
 for i, charmdaughterslist in enumerate(charmdaughters):
  hdecaylen.Fill(decaylen(charmhadrons[i], charmdaughterslist[0]))
  for charmdaughter in charmdaughterslist:
   hkink.Fill(Kinkangle(charmhadrons[i],charmdaughterslist))

 return tracksID

fileinput = r.TFile.Open("ship.conical.Pythia8CharmOnly-TGeant4.root")
inputtree = fileinput.Get("cbmsim")

nevents = 1
#Loop on the events
for ievent in range(1000):
 if (ievent%100==0):
  print ievent
 getdaughtertracks(inputtree,ievent)

#drawing histograms
c0 = r.TCanvas()
hipcharm.Draw()

c1 = r.TCanvas()
hdecaylen.Draw()

c2 = r.TCanvas()
hkink.Draw()
