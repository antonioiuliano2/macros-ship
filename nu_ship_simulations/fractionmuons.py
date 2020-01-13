import ROOT as r
import sys


#how many muons arrive in each detector?
simfile = r.TFile.Open(sys.argv[1])
simtree = simfile.Get("cbmsim")

hrpcid = r.TH1I("hrpcid","Number of hits by initial muon per stations;istation",12,1,12)

hppt = r.TH2D("hppt","P vs PT of initial neutrino;Pt[GeV/c];P[GeV/c]",100,0,10,400,0,400)
hthetap = r.TH2D("hthetap","Angle vs P of initial neutrino;P[GeV/c];#Theta[rad]",400,0,400,100,0,1) 
#starting from target
targetdx = 40.
targetdy = 40.

intarget = 0
print ("{} neutrinos simulated".format(simtree.GetEntries()))
for event in simtree:
 tracks = event.MCTrack
 startneutrino = tracks[0]
 vx = startneutrino.GetStartX() #startx of interaction (the event, i.e. neutrino interaction, starts here)
 vy = startneutrino.GetStartY()
 #angle/momentum
 p = startneutrino.GetP()
 pt = startneutrino.GetPt()
 theta = r.TMath.ASin(pt/p)
 
 hppt.Fill(pt,p)
 hthetap.Fill(p,theta)
 if (r.TMath.Abs(vx) > 40.) or  (r.TMath.Abs(vy) > 40.): #out of target, go to next
  continue 

 intarget = intarget + 1

 rpcpoints = event.ShipRpcPoint

 for point in rpcpoints:
  detID = point.GetDetectorID()
  pdgcode = point.PdgCode()
  trackID = point.GetTrackID()

  #following initial muon
  if (trackID == 1) and (pdgcode == 13):
   hrpcid.Fill(detID-10000)


print("{} neutrinos interacting within target".format(intarget))
#drawing histograms

crpc = r.TCanvas()
hrpcid.Draw()

#canvas out of the way
r.gStyle.SetStatX(0.5);
r.gStyle.SetStatY(0.7);

cangles = r.TCanvas()
cangles.Divide(1,2)
cangles.cd(1)
hppt.Draw("COLZ")
cangles.cd(2)
hthetap.Draw("COLZ")
