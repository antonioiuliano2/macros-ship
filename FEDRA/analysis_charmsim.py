import ROOT 
import fedrarootlogon
import GiuliOpera
#from fedrautils import buildtracks
import recognizecharmdaughters
from collections import Counter #to find most common iterations in a list
from collections import defaultdict # to create a dictionary without knowing the keys in advance

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
hipcharm = ROOT.TH1D("hipcharm","Impact parameter charm daughters",1000,0,1000);
hip = ROOT.TH1D("hip","Impact parameter",1000,0,1000);
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

def createdictionary(tracks):
 d = defaultdict(list)
 for fedratrack in tracks:
  mcevent = fedratrack.MCEvt()
  d[mcevent].append(fedratrack) 

 return d

def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def isinvertex(fedratrack,fedravertex):
  check = False
  trackid = fedratrack.GetSegmentFirst().Track()
  for i in range (fedravertex.N()):
   vertextrack = fedravertex.GetTrack(i)
   vertextrackid = vertextrack.GetSegmentFirst().Track()
   if trackid == vertextrackid:
     check = True
  return check

vertexfile = ROOT.TFile.Open("/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/b000001/17_05_19/vertices_MC_small_modified.root")
vertexrec = vertexfile.Get("EdbVertexRec")

simfile = ROOT.TFile.Open("/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/ship.conical.Pythia8CharmOnly-TGeant4.root")
simtree = simfile.Get("cbmsim")

vertexlist = vertexrec.eVTX
vertextree = vertexfile.Get("vtx")
#all the tracks read ONCE, then create a dictionary
fedratracks = buildtracks("/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/b000001/b000001.0.1.0.trk.root")

tracksperevent = createdictionary(fedratracks)
#loop on the distribution tree
'''
tracks39 = tracksperevent[39]

for t in tracks39:
 print t.eX," ",t.eY," ",t.eZ," ",t.MCEvt()
'''

charminprimary = 0
for event in vertextree:

 ntracks = event.n
 vID = event.vID

 mostprobableMCEventID = Most_Common(event.MCEventID)
 mostprobablemotherid = Most_Common(event.MCMotherID)

 fedravertex = vertexlist.At(vID)
 
 
 charmdaughters = recognizecharmdaughters.getdaughtertracks(simtree,mostprobableMCEventID)
 if mostprobablemotherid == -1:
  #loop on all tracks from the same MCEvent of the vertex
  for track in tracksperevent[mostprobableMCEventID]:

   fedraimpactparameter =  fedravertex.CheckImp(track) #DOES NOT WORK, due to libvt+::VERTEX objects not been saved
   
   elenaimpactparameter = GiuliOpera.IPtoVertex(fedravertex,track)
#   print "Comparing results: ", fedraimpactparameter, elenaimpactparameter
   if track.MCTrack() in charmdaughters:

    hipcharm.Fill(elenaimpactparameter)
    if (isinvertex(track,fedravertex)):
     charminprimary += 1
   else: 
    hip.Fill(elenaimpactparameter)

c0 = ROOT.TCanvas()
hipcharm.Draw()
c1 = ROOT.TCanvas()
hip.Draw()

print "Su un numero di figlie di charm tracciate ",hipcharm.Integral(), " sono state associate al vertice primario ", charminprimary
