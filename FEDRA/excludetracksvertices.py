import ROOT as r
import fedrarootlogon
import fedrautils

trackfilepath = "/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/b000001/b000001.0.1.0.trk.root" #file with all the tracks
vertexfilepath = "/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/b000001/17_05_19/vertices_MC_small_modified.root" #file used for the vertices


dproc = r.EdbDataProc()
gAli = dproc.PVR()
tracklist = fedrautils.buildtracks(trackfilepath, dproc, gAli)

vertexfile = r.TFile.Open(vertexfilepath,'read')
vertexrec = vertexfile.Get("EdbVertexRec")
vertextree = vertexfile.Get("vtx")

vertexlist = vertexrec.eVTX

vertex = vertexlist[0]

ntracks = vertex.N()

for vertex in vertexlist: #loop on vertices
 for itrk in range(ntracks): #loop on tracks from each vertex
  vertextrack = vertex.GetTrack(itrk)
  trackID = vertextrack.GetSegmentFirst().Track()

  tracklist.RemoveAt(trackID)


xv = 0.
yv = 0.

newtracklist = r.TObjArray(100000)

for track in tracklist:
 if (track):
  newtracklist.Add(track)


dproc.MakeTracksTree(newtracklist, xv,yv,"verticesandtracks.root")

#prepare a file with all information for Valerio

newfile = r.TFile.Open("verticesandtracks.root","UPDATE")

newvertexrec = r.EdbVertexRec(vertexrec)

newvertexrec.Write()

#simfile = r.TFile.Open("/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH6_charm/ship.conical.Pythia8CharmOnly-TGeant4.root","read")
#simtree = simfile.Get("cbmsim")

newvertextree = vertextree.CloneTree()
newvertextree.Write()
newsimtree.Write()
newfile.Close()
