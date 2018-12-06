import ROOT as r
import rootlogon
c = r.TCanvas()
def standard_read():
	#opening the file
	f = r.TFile.Open("linked_tracks.root")
	tree = f.Get('tracks')

	tree.SetAlias("t","t.") #setting the alias for the track branch

	#for i in range(tree.GetEntries()): 
	tree.GetEntry(0)
	track = tree.t
	segments = tree.s

	print ('Track TX angle:', track.TX(), "Segment 0 position", segments[0].X(), segments[0].Y(), segments[0].Z())


df = r.RDataFrame('tracks','linked_tracks.root')
histo = df.Histo1D(('hy','Y positions', 100, 0, 100000), 't.eY')
histo.Draw()

df1 = df.Define("t","t.")


