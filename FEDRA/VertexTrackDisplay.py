#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
from argparse import ArgumentParser

fedratracklists = [10,20,30]
trackcolors = [ROOT.kRed, ROOT.kMagenta, ROOT.kBlue] #so we can set different colors for different tracks
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()

parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="fedrafilename", help="file with fedra tracks and vertices",
                    required=True)
parser.add_argument("-n", "--nevent", dest="eventnumber", help="number of event to display", required=True)

options = parser.parse_args()
vertexnumber = options.eventnumber
fedrafilename = options.fedrafilename

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

tracks = buildtracks(fedrafilename, eventnumber)

vertexfile = ROOT.File.Open(fedrafilename)
vertexrec = vertexfile.Get("EdbVertexRec")
vertexlist = vertexrec.eVTX

drawnvertices = ROOT.TObjectArray(100)
drawntracksfromvertex = ROOT.TObjectArray(10000)

vertex = vertexlist[vertexnumber]
ntracksfromvertex = vertex.N()
#adding tracks and vertices to list to be drawn (only one vertex in this case)

drawnvertices.Add(vertex)
for i in range(ntracksfromvertex):
 vertextrack = vertex.GetTrack(i)
 drawntracksfromvertex.Add(track)

def drawtracks(tracks, charmdaughters):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 ds.SetArrTr( tracks )
 ds.SetArrV(varray)
 #print "{} tracks to display\n".format(tracks.GetEntries() )
 ds.Draw();
 #loop on tracks, find charm daughters and replot them with a different color
 for track in tracks:
  if track.GetSegmentFirst().Track() in fedratrackslist: #note, we need to pass to the segments because track() may be confused with the index of the track in the vertex
   nfoundtrack = fedratrackslist.index(track.GetSegmentFirst().Track())
   ds.TrackDraw(track, trackcolors[nfoundtrack])

ROOT.gStyle.SetPalette(1);
dsname="Charm simulation FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)
drawtracks(tracks,charmdaughters)
