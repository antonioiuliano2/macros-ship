import ROOT as r
import fedrarootlogon
import fedrautils
import sys
from argparse import ArgumentParser #not present in good old nusrv9, but the commands should work in a

parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-new", action='store_true') #for new file format

options = parser.parse_args()

dproc = r.EdbDataProc()
gAli = dproc.PVR()
tracklist = fedrautils.buildtracks(options.tracksfilename, dproc, gAli)

if (options.new): #new format, vertex information saved in tree
 r.gROOT.ProcessLine(".L VertexIO.C")
 r.VertexIO.ReadVertexTree(gAli,options.vertexfilename,"1")
 vertexlist = gAli.eVTX

else:
 vertexfile = r.TFile.Open(vertexfilepath,'read')
 vertexrec = vertexfile.Get("EdbVertexRec")
 vertextree = vertexfile.Get("vtx") 

 vertexlist = vertexrec.eVTX


for vertex in vertexlist: #loop on vertices
 ntracks = vertex.N()
 for itrk in range(ntracks): #loop on tracks from each vertex
  vertextrack = vertex.GetTrack(itrk)
  trackID = vertextrack.Track() #taking from the segment, so it tell us the position in the tracklist!
  zpos = vertex.GetVTa(itrk).Zpos() # 1 if track starts from vertex, 0 if track ends in vertex
  if zpos == 1:
    tracklist.RemoveAt(trackID)


xv = 0.
yv = 0.

newtracklist = r.TObjArray(100000)

for track in tracklist:
 if (track):
  newtracklist.Add(track)


dproc.MakeTracksTree(newtracklist, xv,yv,"remainingtracks.root")

#prepare a file with all information for Valerio

#newfile = r.TFile.Open("verticesandtracks.root","UPDATE")

#newvertextree = vertextree.CloneTree() #copying only the tree now, avoiding handling FEDRA objects too much
#newvertextree.Write()
#print "Finished Copying the objects"
#newfile.Close()
