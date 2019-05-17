import ROOT 
import fedrarootlogon
#from fedrautils import buildtracks
import recognizecharmdaughters
from collections import Counter #to find most common iterations in a list
from collections import defaultdict # to create a dictionary without knowing the keys in advance

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
hipcharm = ROOT.TH1D("hipcharm","Impact parameter charm daughters",10000,0,10000);
hip = ROOT.TH1D("hip","Impact parameter",10000,0,10000);
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
for event in vertextree:

 ntracks = event.n
 vID = event.vID

 mostprobableMCEventID = Most_Common(event.MCEventID)

 fedravertex = vertexlist.At(vID)
 
 
 charmdaughters = recognizecharmdaughters.getdaughtertracks(simtree,mostprobableMCEventID)
 if event.n > 8:
  #loop on all tracks from the same MCEvent of the vertex
  for track in tracksperevent[mostprobableMCEventID]:

   impactparameter =  fedravertex.CheckImp(track)
   if track.MCTrack() in charmdaughters:
    hipcharm.Fill(impactparameter)
   else: 
    hip.Fill(impactparameter)

c0 = ROOT.TCanvas()
hipcharm.Draw()
c1 = ROOT.TCanvas()
hip.Draw()

  

