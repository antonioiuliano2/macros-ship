import ROOT as r
import numpy as np

nbinsx = 15
nbinsy = 10
x_bin = {}
q2_bin = {}
for ibin in range((nbinsx*nbinsy)+1):
 x_bin[ibin+1] = []
 q2_bin[ibin+1] = []

expectednumubar_yield = 5.8e+5
eventlist = open("ineccgenieevtnumber_numubarccdis.txt","r")

geniefile = r.TFile.Open("/eos/user/a/aiuliano/public/sims_FairShip/GenieEvents_SHIP/GenieEventsTungsten_01_2023/CCDIS/nu_mu_bar/genie-nu_mu_bar.root")
genietree = geniefile.Get("gst")

hq2_x = r.TH2D("hq2_x",";log10(x);log10(Q2)",nbinsx,-3,0,nbinsy,-1.5,3.5)
#startingloop
for line in eventlist:
 ievent = int(line)
 genietree.GetEntry(ievent)
 logx = r.TMath.Log10(genietree.x)
 logq2 = r.TMath.Log10(genietree.Q2)
 #fill histogram
 hq2_x.Fill(logx, logq2)
 #array of values within bin
 ibin = hq2_x.FindBin(logx, logq2) 
 if(ibin==0):
  print("ibin 0 for values ", logx, "logq2")
  continue
 x_bin[ibin].append(logx)
 q2_bin[ibin].append(logq2)


#draw histogram
r.gStyle.SetOptStat(0)
cq2_x = r.TCanvas()
hq2_x.Scale(1./hq2_x.Integral() * expectednumubar_yield)
hq2_x.Draw("COLZ")
cq2_x.SetLogz()
#cq2_x.SetLogy()
 
gmean = r.TGraph()
 
#mean of arrays
for ibin in range((nbinsx*nbinsy)+1):
 xmean = np.mean(np.array(x_bin[ibin+1]))
 q2mean = np.mean(np.array(q2_bin[ibin+1]))

 print(ibin,xmean,q2mean)
 if(not np.isnan(xmean)):
  gmean.AddPoint(xmean,q2mean)

gmean.SetMarkerColor(r.kRed)
gmean.SetMarkerStyle(r.kFullCircle)
gmean.Draw("P")

