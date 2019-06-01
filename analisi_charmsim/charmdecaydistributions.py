#script to recognize daughters of charm and other particles from the primary proton interaction (created on 13 May 2019)

import ROOT as r
from rootUtils import bookHist
import sys


#opening the file and getting the tree
histofile = r.TFile("distributions_mctrue.root","RECREATE")
hipcharm = r.TH1D("hipcharm","Impact parameter charm daughters",100,0,1000)
hdecaylen= r.TH1D("hdecaylen","Decay length charmed hadrons",1000,0,10000)
hxytovertex = r.TH2D("hxytovertex", "Distance between charmed daughters and vertex", 2000, -1000, 1000, 2000, -1000, 1000)
hztovertex = r.TH1D("hztovertex","Distance between z of charmed daughters and vertex",200,0,20000)
hkink= r.TH1D("hkink","Kink angle charm and daughters",100,0,1)
hkinkprong = []
maxnprongs = 7

#TH1D for recognizing plate position
zstart =  121.8880
zend =  125.5620
hplatez = r.TH1D("hplatez","z Positions, bins as different plates",28,  zstart,  zend)

#for recognizing particles name
pdgdatabase = r.TDatabasePDG.Instance()
#we want to save histograms for the different pdgs
commoncharm = [421, 411, 431, 4122]
hzd0 = r.TH1D("hzd0","Distance between z of charmed daughters and vertex for D0",100,0,100000);
hzdcharge = r.TH1D("hzdcharge","Distance between z of charmed daughters and vertex for D+",100,0,100000);
hzdscharge = r.TH1D("hzdscharge","Distance between z of charmed daughters and vertex for Ds+",100,0,100000);
hzlambdac = r.TH1D("hzlambda+","Distance between z of charmed daughters and vertex for lambac++",100,0,100000);

hpd0 = r.TH1D("hpd0","Momentum of D0",40,0,400);
hpdcharge = r.TH1D("hpdcharge","Momentum of D+",40,0,400);
hpdscharge = r.TH1D("hpdscharge","Momentum of Ds+",40,0,400);
hplambdac = r.TH1D("hplambda+","Momentum of Dlambdac++",40,0,400);

hgammad0 = r.TH1D("hgammad0","Gamma of D0",20,0,200);
hgammadcharge = r.TH1D("hgammadcharge","Gamma of D+",20,0,200);
hgammadscharge = r.TH1D("hgammadscharge","Gamma of Ds+",20,0,200);
hgammalambdac = r.TH1D("hgammalambda+","Gamma of Dlambdac++",20,0,200);

hzcharm = {421: hzd0, 411: hzdcharge, 431: hzdscharge, 4122:hzlambdac}
hpcharm = {421: hpd0, 411: hpdcharge, 431: hpdscharge, 4122:hplambdac}
hgammacharm = {421: hgammad0, 411: hgammadcharge, 431: hgammadscharge, 4122:hgammalambdac}

charmlongntuple = r.TNtuple("charmdecays","Charm Decays","momentum:gamma:pdg:dx:dy:dz:longdecay")

def findplateID(zposition):
	iplate = hplatez.FindBin(zposition)
	return iplate

for iprong in range(maxnprongs):
 hkinkprong.append(r.TH1D("hkink{}prong".format(iprong+1),"Kink angle for {}-prong decay".format(iprong+1), 100, 0, 1))

def decaylen(parenttrack, daughtertrack):

 daughterstartx = daughtertrack.GetStartX()
 daughterstarty = daughtertrack.GetStartZ()
 daughterstartz = daughtertrack.GetStartY()


 parentstartx = parenttrack.GetStartX()
 parentstarty = parenttrack.GetStartZ()
 parentstartz = parenttrack.GetStartY()

 flightlength= pow(daughterstartx-parentstartx,2)+ pow(daughterstarty-parentstarty,2)+ pow(daughterstartz-parentstartz,2)
 
 return pow(flightlength,0.5)
def Kinkangle(parenttrack, daughtertrackslist):

 kink=0.

 parenttx = parenttrack.GetPx()/parenttrack.GetPz()
 parentty = parenttrack.GetPy()/parenttrack.GetPz()
 #computing the average kink angle between mother and track particles
 for daughtertrack in daughtertrackslist:

  daughtertx = daughtertrack.GetPx()/daughtertrack.GetPz()
  daughterty = daughtertrack.GetPy()/daughtertrack.GetPz()
  kink += r.TMath.ATan(r.TMath.Sqrt(pow(parenttx - daughtertx,2) + pow(parentty - daughterty,2)))
 
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

 #delta = r.TMath.Abs(delta)
 #print "DELTA used for impact parameter", delta
 ip = 0.
 vt = 0.
 for i in range(3):
  vt += pow(vertexpos(i) -trackstartpos(i),2) #trying with Pythagora for an independent check
  #ip += pow((vertexpos(i) - trackstartpos(i) - (delta * trackdir(i)/momentum)),2)
  #print "check IP computation", vertexpos(i) - trackstartpos(i), ip
 vt = r.TMath.Sqrt(vt)
 ip = pow(pow(vt,2) - pow(delta,2),0.5)
 #print "FINAL CHECK: ", vt, delta, ip
 return r.TMath.Sqrt(ip)


def getdaughtertracks(inputtree,eventnumber):

 #interesting pdgcodes to check
 intermediatelist = [223, 3332, 3224, 331, 221, 20213, 3212, 213, 113, 2224, 323] #particles with very short lifetime
 signallist = [431, 411, 4122, 421, 4132, 4232, 4332, 441] #charmed hadrons

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

   if pdgdatabase.GetParticle(pdgcode):
    name = pdgdatabase.GetParticle(pdgcode).GetName()
    charge = pdgdatabase.GetParticle(pdgcode).Charge()

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
     if momentum > 0.1 and r.TMath.Abs(charge) > 0: #energy cut, also we select charged particles
  
      charmindex = charmpdgs.index(motherpdg) #what charm was found?
      charmdaughters[charmindex].append(track)

      #print "Track Momentum", track.GetPx(), track.GetPy(), track.GetPz()
      #print "Start of Track: ",track.GetStartX(), track.GetStartY(), track.GetStartZ()
      #print "Vertex position: ",vertex(0), vertex(1), vertex(2)
      impactparameter = IPtoVertex(vertex,track)
      #print "Estimated IP ", impactparameter, "in micron ", impactparameter*cmtomicron
      hipcharm.Fill(impactparameter*cmtomicron)

      tracksID.append(j)

 #loop on the two groups of charmdaughters
 for i, charmdaughterslist in enumerate(charmdaughters):
  if len(charmdaughterslist) > 0: #we need at least one charged daughter to fill the histograms
   #for charmdaughter in charmdaughterslist:
   averagekinkangle = Kinkangle(charmhadrons[i],charmdaughterslist)
   hkink.Fill(averagekinkangle)
   hxytovertex.Fill((charmdaughterslist[0].GetStartX() - vertex(0))*cmtomicron,(charmdaughterslist[0].GetStartY() - vertex(1))*cmtomicron)
   hztovertex.Fill((charmdaughterslist[0].GetStartZ() - vertex(2))*cmtomicron)
   hdecaylen.Fill(decaylen(charmhadrons[i],charmdaughterslist[0])*cmtomicron)	
   charmhadronpdg = charmhadrons[i].GetPdgCode()
   mass = pdgdatabase.GetParticle(charmhadronpdg).Mass()
   gamma = charmhadrons[i].GetEnergy()/mass
        
         
   longdecay = True
   primaryvertexplate = findplateID(vertex(2))
   secondaryvertexplate = findplateID(charmdaughterslist[0].GetStartZ())
   if(primaryvertexplate==secondaryvertexplate): longdecay = False
   dx = (charmdaughterslist[0].GetStartX() - vertex(0))*cmtomicron
   dy = (charmdaughterslist[0].GetStartY() - vertex(1))*cmtomicron
   dz = (charmdaughterslist[0].GetStartZ() - vertex(2))*cmtomicron
   momentum = charmhadrons[i].GetP()
   pdgcode = charmhadronpdg
         
   charmlongntuple.Fill(momentum,gamma,pdgcode,dx,dy,dz,longdecay)

   if r.TMath.Abs(charmhadronpdg) in commoncharm:
         hzcharm[r.TMath.Abs(charmhadronpdg)].Fill((charmdaughterslist[0].GetStartZ() - vertex(2))*cmtomicron)
         hpcharm[r.TMath.Abs(charmhadronpdg)].Fill(charmhadrons[i].GetP())
         hgammacharm[r.TMath.Abs(charmhadronpdg)].Fill(gamma)
   #for charmdaughter in charmdaughterslist:    
         

   #print "TEST ",len(charmdaughterslist)
   hkinkprong[len(charmdaughterslist)-1].Fill(averagekinkangle)
 return tracksID

fileinput = r.TFile.Open("ship.conical.Pythia8CharmOnly-TGeant4.root")
inputtree = fileinput.Get("cbmsim")

nevents = 1
#Loop on the events
for ievent in range(1000):
 print "Start of event: ", ievent
 getdaughtertracks(inputtree,ievent)

#drawing histograms and saving them to file
histofile.cd("")
c0 = r.TCanvas()
hipcharm.GetXaxis().SetTitle("ip[#mum]")
hipcharm.Draw()
hipcharm.Write()

c1 = r.TCanvas()
hdecaylen.GetXaxis().SetTitle("decaylength[cm]")
hdecaylen.Draw()
hdecaylen.Write()

cdistance = r.TCanvas()
cdistance.Divide(1,2)
cdistance.cd(1)
hxytovertex.GetXaxis().SetTitle("x[#mum]")
hxytovertex.GetYaxis().SetTitle("y[#mum]")
hxytovertex.Draw("COLZ")
hxytovertex.Write()
cdistance.cd(2)
hztovertex.Draw()
hztovertex.GetXaxis().SetTitle("z[#mum]")
hztovertex.Write()

czcharm = r.TCanvas()
czcharm.Divide(2,2)
czcharm.cd(1)
hzcharm[421].Draw()
hzcharm[421].GetXaxis().SetTitle("z[#mum]")
hzcharm[421].Write()
czcharm.cd(2)
hzcharm[411].Draw()
hzcharm[411].GetXaxis().SetTitle("z[#mum]")
hzcharm[411].Write()
czcharm.cd(3)
hzcharm[431].Draw()
hzcharm[431].GetXaxis().SetTitle("z[#mum]")
hzcharm[431].Write()
czcharm.cd(4)
hzcharm[4122].Draw()
hzcharm[4122].GetXaxis().SetTitle("z[#mum]")
hzcharm[4122].Write()
czcharm.Write()

cpcharm = r.TCanvas()
cpcharm.Divide(2,2)
cpcharm.cd(1)
hpcharm[421].Draw()
hpcharm[421].GetXaxis().SetTitle("p[GeV]")
hpcharm[421].Write()
cpcharm.cd(2)
hpcharm[411].Draw()
hpcharm[411].GetXaxis().SetTitle("p[GeV]")
hpcharm[411].Write()
cpcharm.cd(3)
hpcharm[431].Draw()
hpcharm[431].GetXaxis().SetTitle("p[GeV]")
hpcharm[431].Write()
cpcharm.cd(4)
hpcharm[4122].Draw()
hpcharm[4122].GetXaxis().SetTitle("p[GeV]")
hpcharm[4122].Write()
cpcharm.Write()

cgammacharm = r.TCanvas()
cgammacharm.Divide(2,2)
cgammacharm.cd(1)
hgammacharm[421].Draw()
hgammacharm[421].GetXaxis().SetTitle("gamma")
hgammacharm[421].Write()
cgammacharm.cd(2)
hgammacharm[411].Draw()
hgammacharm[411].GetXaxis().SetTitle("gamma")
hgammacharm[411].Write()
cgammacharm.cd(3)
hgammacharm[431].Draw()
hgammacharm[431].GetXaxis().SetTitle("gamma")
hgammacharm[431].Write()
cgammacharm.cd(4)
hgammacharm[4122].Draw()
hgammacharm[4122].GetXaxis().SetTitle("gamma")
hgammacharm[4122].Write()
cgammacharm.Write()

charmlongntuple.Write()


c2 = r.TCanvas()
hkink.GetXaxis().SetTitle("rad")
hkink.Draw()
hkink.Write()

c3 = r.TCanvas()
c3.Divide(3,3)
for iprong in range(maxnprongs):
 c3.cd(iprong+1)
 hkinkprong[iprong].SetTitle("rad")
 hkinkprong[iprong].Draw()
 hkinkprong[iprong].Write()
c3.Write()
