from __future__ import division
import ROOT as r
import numpy as np

inputfile = r.TFile.Open("ship.conical.Genie-TGeant4.root")
simtree = inputfile.Get("cbmsim")

hnrpcsmuon = r.TH1I("hnrpcsmuon","Number of rpcs transversed by each muon track",11,1,12)
hnrpcsother = r.TH1I("hnrpcsother","Number of rpcs transversed by each not muon track",11,1,12)

nevents = simtree.GetEntries()

def countrpcstations(mytree):
 maxtracksinrpc = int(1E+5)
 nrpcbytrack = np.zeros((maxtracksinrpc,2), dtype=int) #how many rpc points for each track in this event
 
 rpcpoints = mytree.ShipRpcPoint

 for point in rpcpoints:
  trackID = point.GetTrackID()
  if trackID > 0: 
   nrpcbytrack[trackID,1] = nrpcbytrack[trackID,1] +1

 return nrpcbytrack

for ievent in range(nevents):
 if ievent % 100 == 0:
  print "Arrivato all'evento ", ievent
 simtree.GetEntry(ievent)
 tracks = simtree.MCTrack
 reftrack = tracks[0]

 inacceptance = False
 if abs(reftrack.GetStartX()) > 40. or abs(reftrack.GetStartY()) > 40.:
  inacceptance = True
 weight = reftrack.GetWeight()

 nrpcbytrack = countrpcstations(simtree)

 for irow in range(nrpcbytrack.shape[0]):
  if nrpcbytrack[irow,1]>0: 
   if (abs(tracks[irow].GetPdgCode())==13):
    hnrpcsmuon.Fill(nrpcbytrack[irow,1])
   else:
    hnrpcsother.Fill(nrpcbytrack[irow,1])

hnrpcsother.Draw()
hnrpcsmuon.SetLineColor(r.kRed)
hnrpcsmuon.Draw("SAMES")

def muonprobability(hnrpcsother, hrpcsmuon):

 nminrpcs = 7; #cut for muon identification
 nrpcsmuons = hnrpcsmuon.Integral() #how many muons passing in first rpc
 nrpcsother = hnrpcsother.Integral() #how many other particles passing in first rpc
 fracmuons = nrpcsmuons/(nrpcsmuons+nrpcsother)
 fracother = nrpcsother/(nrpcsmuons+nrpcsother)

 efficiency = hnrpcsmuon.Integral(nminrpcs,12)/hnrpcsmuon.Integral()
 misidentification = hnrpcsother.Integral(nminrpcs,12)/hnrpcsother.Integral()

 probmuon = efficiency * fracmuons /(efficiency * fracmuons + misidentification * fracother) #applicatioo of Bayes theorem for probability of being a muons, given a positive response

 print "Fraction of muons: {:.2f}".format(fracmuons)
 print "Efficiency of at least {} RPCs stations: {:.2f}, misidentification: {:.4f}".format(nminrpcs, efficiency, misidentification)
 print "Probability of being a muon, if it passes at least {} RPCs stations: {:.4f}".format(nminrpcs, probmuon)

muonprobability(hnrpcsother, hnrpcsmuon)
