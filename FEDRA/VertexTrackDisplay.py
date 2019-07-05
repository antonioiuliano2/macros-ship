#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
#usage: python -i VerteTrackDisplay.py inputfile nevent

#from argparse import ArgumentParser #not present in good old nusrv9, but the commands should work in a reasonable python setup, only need to remove the parser and options comments,then comment the sys.argv lines

fedratrackslist = [10,20,30]
trackcolors = [ROOT.kRed, ROOT.kMagenta, ROOT.kBlue] #so we can set different colors for different tracks
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()

#parser = ArgumentParser()
#parser.add_argument("-f", "--fedra", dest="fedrafilename", help="file with fedra tracks and vertices",
     #               required=True)
#parser.add_argument("-n", "--nevent", dest="eventnumber", help="number of event to display", required=True)

#options = parser.parse_args()
#vertexnumber = options.eventnumber
#fedrafilename = options.fedrafilename

fedrafilename = sys.argv[1]
vertexnumber = int(sys.argv[2])

def buildtracks(filename):
 
 dproc.ReadTracksTree(gAli,filename,"nseg>1")
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

tracks = buildtracks(fedrafilename)

vertexfile = ROOT.TFile.Open(fedrafilename)
vertexrec = vertexfile.Get("EdbVertexRec")
vertexlist = vertexrec.eVTX

drawnvertices = ROOT.TObjArray(100)
drawntracksfromvertex = ROOT.TObjArray(10000)

vertex = vertexlist.At(vertexnumber)
ntracksfromvertex = vertex.N()
#adding tracks and vertices to list to be drawn (only one vertex in this case)

drawnvertices.Add(vertex)
for i in range(ntracksfromvertex):
 vertextrack = vertex.GetTrack(i)
 drawntracksfromvertex.Add(vertextrack)

def drawtracks(vertextracks,tracks):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 ds.SetArrTr( vertextracks )
 ds.SetArrV(drawnvertices)
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
drawtracks(drawntracksfromvertex,tracks)
