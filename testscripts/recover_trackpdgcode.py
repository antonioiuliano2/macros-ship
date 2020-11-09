#need to recover trackpdg
import ROOT as r
import fedrarootlogon

vertexfile = r.TFile.Open("vertextree_test.root")
vertextree = vertexfile.Get("vtx")

trackfile = r.TFile.Open("linked_tracks.root")
tracktree = trackfile.Get("tracks")

tracktree.SetAlias("track","t.")
#start vertex  loop
for ivtx, vtx in enumerate(vertextree):

 trackIDs = vtx.TrackID
 #start track loop
 for tID in trackIDs:
  tracktree.GetEntry(tID)
  recotrack = tracktree.track
  MCPdgCode = recotrack.Vid(0)
 if (ivtx % 1000 == 0):
  print ("Arrived at entry, ",ivtx)
  
