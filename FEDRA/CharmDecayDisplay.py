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
parser.add_argument("-n", "--nevent", dest="eventnumber", help="number of event to display", required=True)

options = parser.parse_args()
eventnumber = options.eventnumber
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
tracks = buildtracks(fedrafilename, eventnumber)


def drawtracks(tracks, charmdaughters):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 ds.SetArrTr( tracks )
 print "{} tracks to display\n".format(tracks.GetEntries() )
 ds.Draw();
 #loop on tracks, find charm daughters and replot them with a different color
 for track in tracks:
  if track.MCTrack() in charmdaughters:
   ds.TrackDraw(track, ROOT.kMagenta)
  if track.MCTrack() in charmIDs:
   ds.TrackDraw(track, ROOT.kBlue)
   testrmax = decaysearch.FedraTrackKink(track)
   print "Prova ", testrmax

ROOT.gSystem.Load("/afs/cern.ch/work/a/aiuliano/public/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so")
decaysearch = ROOT.ShipCharmDecaySearch()
simfile = ROOT.TFile.Open(shipfilename)
cbmsim = simfile.Get("cbmsim")
(charmIDs, charmdaughters) = recognizecharmdaughters.getdaughtertracks(cbmsim,int(eventnumber))
ROOT.gStyle.SetPalette(1);
dsname="Charm simulation FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)
drawtracks(tracks,charmdaughters)
