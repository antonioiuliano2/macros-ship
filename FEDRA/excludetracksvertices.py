import ROOT as r
import fedrarootlogon
import fedrautils

trackfilepath = "linked_tracks.root" #file with all the tracks
vertexfilepath = "vertices_firstquarter.root" #file used for the vertices


dproc = r.EdbDataProc()
gAli = dproc.PVR()
tracklist = fedrautils.buildtracks(trackfilepath, dproc, gAli)

vertexfile = r.TFile.Open(vertexfilepath,'read')
vertexrec = vertexfile.Get("EdbVertexRec")
vertextree = vertexfile.Get("vtx")

vertexlist = vertexrec.eVTX

vertex = vertexlist[0]

for vertex in vertexlist: #loop on vertices
 ntracks = vertex.N()
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
newvertexrec.eVTX = vertexrec.eVTX

newvertexrec.Write()

newvertextree = vertextree.CloneTree()
newvertextree.Write()
print "Finished Copying the objects"
newfile.Close()
