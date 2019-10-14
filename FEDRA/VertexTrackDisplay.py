#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
#usage: python -i VerteTrackDisplay.py -f vertexfile -t tracksfile -nt trackIDS -nv vertexIDS

from argparse import ArgumentParser #not present in good old nusrv9, but the commands should work in a reasonable python setup, only need to remove the parser and options comments,then comment the sys.argv lines

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()

fedratrackslist = []
#vertexnumberlist = [10, 20]
isolatedtrackcolors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta] #so we can set different colors for different tracks
vertextrackcolors = [ROOT.kYellow,ROOT.kBlue,ROOT.kRed] #so we can set different colors for different tracks
#list of possible options
parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-nv", "--nvertices", nargs='+', dest="vertexnumberlist", help="number of vertices to display", required=True)
parser.add_argument("-nt", "--ntracks", nargs='*', dest="tracklist", help="number of isolated tracks to display")
parser.add_argument("-new", action='store_true') #for new file format

options = parser.parse_args()
vertexnumberlist = options.vertexnumberlist
if (options.tracklist): 
 fedratrackslist = options.tracklist
 fedratrackslist = map(int, fedratrackslist)
vertexfilename = options.vertexfilename
tracksfilename = options.tracksfilename


inputfile = ROOT.TFile.Open(tracksfilename,"READ")
tracktree = inputfile.Get("tracks")
tracktree.SetAlias("trk","t.") #points create confusion to python

tracktree.BuildIndex("trid")
tracks = []
for trackID in fedratrackslist:
 tracktree.GetEntryWithIndex(trackID)
 #temporary object for reading the file and building EdbTrackP
 temptrack = ROOT.EdbTrackP()
 temptrack .Copy(ROOT.EdbTrackP(tracktree.trk))
 segments = tracktree.s
 fittedsegments = tracktree.sf
 #loop on segments associated to the track
 for seg in segments:
     temptrack .AddSegment(seg)
     #temptrack .AddSegmentF(segf)
     temptrack .SetSegmentsTrack(temptrack.ID())
     temptrack .SetCounters()
 mytrack = ROOT.EdbTrackP()
 mytrack.Copy(temptrack)
 tracks.append(mytrack)

drawnvertices = ROOT.TObjArray(100)
drawntracksfromvertex = ROOT.TObjArray(10000)

if (options.new): #new format, vertex information saved in tree
 ROOT.gROOT.ProcessLine(".L VertexIO.C")
 for vertexnumber in vertexnumberlist:
  vertex = ROOT.VertexIO.GetVertexFromTree(gAli,vertexfilename,int(vertexnumber))
  ntracksfromvertex = vertex.N()
  #adding tracks and vertices to list to be drawn (only one vertex in this case)
  drawnvertices.Add(vertex)

  for i in range(ntracksfromvertex):
   vertextrack = vertex.GetTrack(i)
   drawntracksfromvertex.Add(vertextrack)
   #dproc.ReadVertexTree(gAli,vertexfilename,"nseg>1")

else:
 vertexfile = ROOT.TFile.Open(vertexfilename)
 vertexrec = vertexfile.Get("EdbVertexRec")
 vertexlist = vertexrec.eVTX

#adding vertices
 for vertexnumber in vertexnumberlist:
  vertex = vertexlist.At(int(vertexnumber))
  ntracksfromvertex = vertex.N()
#adding tracks and vertices to list to be drawn (only one vertex in this case)
  drawnvertices.Add(vertex)

  for i in range(ntracksfromvertex):
   vertextrack = vertex.GetTrack(i)
   drawntracksfromvertex.Add(vertextrack)

def drawtracks(vertextracks,othertracks):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 ds.SetArrTr( vertextracks )
 ds.SetArrV(drawnvertices)
 ds.Draw()
 print len(othertracks),"other tracks to display\n"
 for itrk, track in enumerate(othertracks):
   ds.TrackDraw(track,isolatedtrackcolors[itrk])
 #loop on vertices to draw associated tracks
 for ivtx, vertex in enumerate(drawnvertices):
  for itrk in range(vertex.N()):#tracks associated to that vertex
   track = vertex.GetTrack(itrk)
   ds.TrackDraw(track, vertextrackcolors[ivtx])

ROOT.gStyle.SetPalette(1);
dsname="Charm simulation FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)
drawtracks(drawntracksfromvertex,tracks)

