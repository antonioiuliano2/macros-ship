import ROOT as r
import fedrarootlogon

dproc = r.EdbDataProc()
gAli = dproc.PVR()

def buildtracks(filename,cut):
 print cut
 dproc.ReadTracksTree(gAli,filename,cut)
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks
 
 return tracks

xmin = 20000.
ymin = 20000.
xmax = 90000.
ymax = 90000.

nnumbers = 100

def directextraction(nnumbers):
 for i in range(nnumbers):
  #generating a position in the space
  vx = r.gRandom.Uniform(xmin,xmax)
  vy = r.gRandom.Uniform(ymin,ymax)

  print "generated point in ", vx, vy
  #reading tracks in a 6 mm x 6 mm region around this point
  tracklist = buildtracks("/home/antonio/Dottorato/Charmdata/CH1-R6/charm1_run6_opera1_complete.trk.root",("nseg>1&&TMath::Abs(s[0].eX-%d)<6000&&TMath::Abs(s[0].eY-%d)<6000"%(vx,vy)))

  gAli.ResetTracks()

def createdictionary(tracks,hreticule):
 d = defaultdict(list)
 for fedratrack in tracks:
  #taking the start coordiantes
  startx = fedratrack.GetSegmentFirst().X()
  starty = fedratrack.GetSegmentFirst().Y()

  binnumber = hreticule.FindBin(startx,starty)

  d[binnumber].append(fedratrack) 

 return d

def reticule(nnumbers):
 tracklist = buildtracks("/home/antonio/Dottorato/Charmdata/CH1-R6/charm1_run6_opera1_complete.trk.root","nseg>1")
 
 hreticule = ROOT.TH2F("hreticule","Map to associate the track dictionary",125,0,125000,100,0,10000)
 tracksdictionary = createdictionary(tracklist,hreticule)

 #starting loop
 for i in range(nnumbers):
   #generating a position in the space
   vx = r.gRandom.Uniform(xmin,xmax)
   vy = r.gRandom.Uniform(ymin,ymax)

   binnumber = hreticule.FindBin(vx,vy)
   tracklist = tracksdictionary[binnumber]

reticule(nnumbers)