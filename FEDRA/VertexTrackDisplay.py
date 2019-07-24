#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys

ROOT.gSystem.Load("/afs/cern.ch/work/a/aiuliano/public/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so")
#usage: python -i VerteTrackDisplay.py inputfile nevent

from argparse import ArgumentParser #not present in good old nusrv9, but the commands should work in a reasonable python setup, only need to remove the parser and options comments,then comment the sys.argv lines

fedratrackslist = []
#vertexnumberlist = [10, 20]
isolatedtrackcolors = [ROOT.kRed, ROOT.kMagenta, ROOT.kYellow, ROOT.kBlue] #so we can set different colors for different tracks
vertextrackcolors = [ROOT.kMagenta,ROOT.kBlue,ROOT.kYellow,ROOT.kOrange,ROOT.kRed] #so we can set different colors for different tracks
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()

parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="fedrafilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-nv", "--nvertices", nargs='+', dest="vertexnumberlist", help="number of vertices to display", required=True)
parser.add_argument("-nt", "--ntracks", nargs='*', dest="tracklist", help="number of isolated tracks to display")

options = parser.parse_args()
vertexnumberlist = options.vertexnumberlist
if (options.tracklist): 
 fedratrackslist = options.tracklist
 fedratrackslist = map(int, fedratrackslist)
fedrafilename = options.fedrafilename
tracksfilename = options.tracksfilename

#fedrafilename = sys.argv[1]
#vertexnumber = int(sys.argv[2])

def buildtracks(filename):
 
 dproc.ReadTracksTree(gAli,tracksfilename,"nseg>1")
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

#adding vertices
for vertexnumber in vertexnumberlist:
 vertex = vertexlist.At(int(vertexnumber))
 ntracksfromvertex = vertex.N()
#adding tracks and vertices to list to be drawn (only one vertex in this case)
 drawnvertices.Add(vertex)

 for i in range(ntracksfromvertex):
  vertextrack = vertex.GetTrack(i)
  drawntracksfromvertex.Add(vertextrack)

graphip = ROOT.TGraphErrors()

def fillip(vertexpos, track):
 decaysearch = ROOT.ShipCharmDecaySearch()
 nseg = track.N()
 ipoint = 0
 print "Test ", vertexpos[0], vertexpos[1], vertexpos[2]
 for iseg in range(nseg):
  segment = track.GetSegment(iseg)
  segpos = ROOT.TVector3(segment.X(),segment.Y(),segment.Z())
  print "ReTest ", segpos[0], segpos[1], segpos[2]
  dz = vertexpos[2] - segpos[2];
  st = 0.003
  sx = 1
  sy = sx
  varipx = dz*dz*st*st+sx*sx
  varipy = dz*dz*st*st+sy*sy
  sigmaip = ROOT.TMath.Sqrt(varipx+varipy)
  graphip.SetPoint(ipoint, segment.Plate(), decaysearch.IPtoVertex(vertexpos,segpos,segment.TX(), segment.TY()))
  graphip.SetPointError(ipoint, 0, sigmaip)
  ipoint = ipoint+1

def drawtracks(vertextracks,tracks):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 ds.SetArrTr( vertextracks )
 ds.SetArrV(drawnvertices)
 ds.Draw()
 #print "{} tracks to display\n".format(tracks.GetEntries() )
 #loop on tracks, find charm daughters and replot them with a different color
 #print "PROVA:", len(tracks)
 for track in tracks:
  if track.GetSegmentFirst().Track() in fedratrackslist: #note, we need to pass to the segments because track() may be confused with the index of the track in the vertex    
   nfoundtrack = fedratrackslist.index(track.GetSegmentFirst().Track())
   ds.TrackDraw(track, isolatedtrackcolors[nfoundtrack])
 #loop on vertices to draw associated tracks
 for ivtx, vertex in enumerate(drawnvertices):
  for itrk in range(vertex.N()):#tracks associated to that vertex
   track = vertex.GetTrack(itrk)
   if track.MCTrack() == 1:
    print "Test "
    vx = vertex.X()
    vy = vertex.Y()
    vz = vertex.Z()
    vertexpos = ROOT.TVector3(vx,vy,vz)
    fillip(vertexpos,track)
    ds.TrackDraw(track,ROOT.kBlue)
   else:
    ds.TrackDraw(track, vertextrackcolors[ivtx])

ROOT.gStyle.SetPalette(1);
dsname="Charm simulation FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)
drawtracks(drawntracksfromvertex,tracks)

c1 = ROOT.TCanvas()
graphip.Draw("AP*")
graphip.GetXaxis().SetTitle("NPlate")
graphip.GetYaxis().SetTitle("ip[#mu m]")

