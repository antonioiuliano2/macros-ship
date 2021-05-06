import ROOT as r
import fedrarootlogon
c = r.TCanvas()

#opening the file
f = r.TFile.Open("linked_tracks.root")
tree = f.Get('tracks')
tree.SetAlias("t","t.") #setting the alias for the track branch
#for i in range(tree.GetEntries()): 
tree.GetEntry(0)
track = tree.t
segments = tree.s
fittedsegments = tree.sf

print ('Track TX angle:', track.TX(), "Segment 0 position", segments[0].X(), segments[0].Y(), segments[0].Z())

for (seg, segf) in zip(segments, fittedsegments):
 #we need to access Covariance matrix
 mycov = segf.COV()
 mycov.Print()


