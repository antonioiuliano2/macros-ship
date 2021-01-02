#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
#usage: python -i VerteTrackDisplay.py -f vertexfile -t tracksfile -nt trackIDS -nv vertexIDS

from argparse import ArgumentParser 

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
gAli.SetScanCond(scancond)

fedratrackslist = []
#vertexnumberlist = [10, 20]
isolatedtrackcolors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta] #so we can set different colors for different tracks
vertextrackcolors = [ROOT.kRed,ROOT.kBlue,ROOT.kRed, ROOT.kGreen, ROOT.kCyan, ROOT.kGray] #so we can set different colors for different tracks
#list of possible options
parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("--shower", action='store_true')
parser.add_argument("-ns","--nshower",dest="nshower",help="number of shower to draw")
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

proc = ROOT.EdbDataProc()

if (options.new): #new format, vertex information saved in tree
 for vertexnumber in vertexnumberlist:

  vertexrec = ROOT.EdbVertexRec()
  vertexrec.SetPVRec(gAli)
  vertexrec.eDZmax=3000.
  vertexrec.eProbMin=0.01
  vertexrec.eImpMax=15.
  vertexrec.eUseMom=False
  vertexrec.eUseSegPar=True
  vertexrec.eQualityMode=0

  vertex = proc.GetVertexFromTree(vertexrec,vertexfilename,int(vertexnumber))
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

if (options.shower): # draw shower
  ROOT.gROOT.ProcessLine(".L /home/antonio/Scrivania/macros-ship/FEDRA/shower_reconstruction/SimpleShowerRecInterface.C")
  showerrecinterface = ROOT.SimpleShowerRecInterface()
  #find which pvrec is it, according to my ordering
  startx = tracks[0].X()
  starty = tracks[0].Y()
  xcode = int(startx/10000) * 10000
  ycode = int(starty/10000) * 10000
  print(xcode, ycode)
  pvrecfile = ROOT.TFile.Open("pvrecs/pvrec_{}_{}.root".format(xcode,ycode))
  showerrecinterface.LoadPVRec(pvrecfile)
  seglist = showerrecinterface.DrawShower(int(options.nshower),"/home/antonio/cernbox/Synched/Charmsim_Showreco/Showreco_ds_25_01/Shower_{}_{}.root".format(xcode,ycode))

def drawtracks(vertextracks,othertracks):
 #ds.SetVerRec(gEVR);
 ds.SetDrawTracks(4)
 if (options.shower):
   #print(seglist[2].X(),seglist[2].Y(),seglist[2].Z(),seglist[2].DZ())
   myseglist = ROOT.TObjArray()
   for segment in seglist:
     mysegline = ds.SegLine(segment)
     mysegline.SetLineColor(ROOT.kBlue)
     mysegline.SetLineWidth(6)
     myseglist.Add(mysegline)
   ds.SetArrSegG(myseglist)
 ds.SetArrTr( vertextracks )
 ds.SetArrV(drawnvertices)
 ds.Draw()
 #loop on vertices to draw associated tracks
 for ivtx, vertex in enumerate(drawnvertices):
  for itrk in range(vertex.N()):#tracks associated to that vertex
   track = vertex.GetTrack(itrk)
   if (ivtx < len(vertextrackcolors)):
     vertexcolor =  vertextrackcolors[ivtx]
   else:
     vertexcolor = ROOT.kWhite
   ds.TrackDraw(track,vertexcolor)
 print (len(othertracks),"other tracks to display\n")
 for itrk, track in enumerate(othertracks):
   ds.TrackDraw(track,isolatedtrackcolors[itrk])
 #if (options.shower):
 #  for segment in seglist:
 #    ds.TrackDraw(segment)   

ROOT.gStyle.SetPalette(1);

dsname="Charm simulation FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)


drawtracks(drawntracksfromvertex,tracks)

#get vertex file and vertex tree
vertexfile = ROOT.TFile.Open(vertexfilename,"read")
vertextree = vertexfile.Get("vtx")

def close():
  '''Close canvas, allowing program to exit without crashing ROOT'''
  global ds
  ds = 0

def drawnearvertices(x,y,z,radius,minmolt):
  '''draw vertices near a given point'''
  nearvertices = vertextree.CopyTree("TMath::Sqrt(pow(vx-{},2)+pow(vy-{},2)+pow(vz-{},2))<{} && n >= {}".format(x,y,z,radius,minmolt))
  print ("Other {} vertices to draw".format(nearvertices.GetEntries()))
  for ivtx,entry in enumerate(nearvertices):
    vertex = proc.GetVertexFromTree(gAli,vertexfilename,entry.vID)
    if str(vertex.ID()) in vertexnumberlist:
       print ("Already present 1 vertex")
       continue
    if (vertex.Flag()==2 or vertex.Flag()==5):
       continue
    ds.VertexDraw(vertex)
    print("Test: ",vertex.VZ(),vertex.Z())
    for itrk in range(vertex.N()):#tracks associated to that vertex
     track = vertex.GetTrack(itrk)
     ds.TrackDraw(track,ivtx+2)
