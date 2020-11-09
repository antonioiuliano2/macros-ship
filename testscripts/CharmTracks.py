#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
import recognizecharmdaughters
from argparse import ArgumentParser

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()

parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="fedrafilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-s", "--ship", dest="shipfilename", help="output from FairShip simulation", required=True)

options = parser.parse_args()
fedrafilename = options.fedrafilename
shipfilename = options.shipfilename

def buildtracks(filename,eventnumber):
 
 dproc.ReadTracksTree(gAli,filename,"nseg>1&&s.eMCEvt=={}".format(eventnumber))
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks

 np = gAli.Npatterns();
 for i in range (np): 
    p = gAli.GetPattern(i);
    ns = p.N();
    for j in range(ns):
      p.GetSegment(j).SetDZ(300);

 return tracks

#start of the main loop script

#inputfilename = "/ship/CHARM2018/CH1-R6/b000001/b000001.0.0.0.trk.root"
simfile = ROOT.TFile.Open(shipfilename)
cbmsim = simfile.Get("cbmsim")
nevents = cbmsim.GetEntries()

trackedcharm = 0
mergedcharm = 0

for ievent in range(nevents):
 (charmIDs, charmdaughters) = recognizecharmdaughters.getdaughtertracks(cbmsim,int(ievent))
 tracks = buildtracks(fedrafilename, ievent)

 for track in tracks:

  if track.MCTrack() in charmIDs:
   mergedtracks = False
   trackedcharm = trackedcharm + 1
   nseg = track.N()
   for iseg in range(nseg):
    segment = track.GetSegment(iseg)
    if (segment.MCTrack() != track.MCTrack()): 
     mergedtracks = True
     print "possible merge: ",segment.MCTrack(), segment.MCEvt(), track.MCTrack(), track.MCEvt()
   if (mergedtracks): mergedcharm = mergedcharm + 1

print "End of loop"
print "Number of charm tracked: {} , from which {} merged with daughters".format(trackedcharm, mergedcharm)
