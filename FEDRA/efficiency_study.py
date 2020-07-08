'''python version of efficiency_studyC, to easily hande multiple files
Usage: 
1)for only one file
  python -i efficiency_study.py linked_tracks.root 
2)for 4 surface quarters 
  python -i efficiency_study.py tracks1.root tracks2.root tracks3.root tracks4.root
created by Antonio on date 16 June 2020
'''
import ROOT as r
import sys
import fedrarootlogon

nplates=57

goodtrack ="t.eFlag>=0 &&t.eProb>0.01&&npl >= 15"

surfaces = []

#need to use splicing to remove default first argument
filenames = sys.argv[1:]
nfiles = len(filenames)
if nfiles==1:
  surfaces.append("1")
elif nfiles==4:
  surfaces.append("t.eX < 62500 && t.eY<50000")
  surfaces.append("t.eX > 62500 && t.eY<50000")
  surfaces.append("t.eX < 62500 && t.eY>50000")
  surfaces.append("t.eX > 62500 && t.eY>50000")
else:
  print("ERROR, please provide 1 or 4 tracks files as input")

#degining histograms
hexpected = r.TH1D("hexpected", "Tracks expected to be found in each plate", nplates,1,nplates+1);
hfound = r.TH1D("hfound", "Tracks with an associated segment in each plate;Plate", nplates,1,nplates+1); 

hxy = r.TH2D("hxy", "2D position of tracks;x[mm];y[mm]",125,0,125,100,0,100);

def trackloop(filename, condition):
 '''loop over all tracks in filename according to condition'''
 #getting file and tree
 firstfile = r.TFile.Open(filename)
 firsttracks = firstfile.Get("tracks")

 #reading entries from selection
 print("Selection of good tracks:")
 print(condition)
 firsttracks.Draw(">>lst",r.TCut(condition))

 lst = r.gDirectory.GetList().FindObject("lst");

 nlst = lst.GetN()

 #doing loop
 print("doing loop over ",nlst," tracks over ", firsttracks.GetEntries())

 for entry in range(nlst):
  ientry = lst.GetEntry(entry)
  firsttracks.GetEntry(ientry)

  #filling histogram with xy start coordinates
  hxy.Fill(firsttracks.s[0].X()/1000., firsttracks.s[0].Y()/1000.) 

  #getting first segment and last segment of track
  nseg = firsttracks.nseg

  firstplate = firsttracks.s[0].Plate()
  lastplate = firsttracks.s[nseg-1].Plate()
  #I expect to found the track in all segments from first to last
  for iplate in range(firstplate, lastplate+1):
   hexpected.Fill(iplate)
  #hexpected.Fill(0)
  #next, found plates
  for seg in firsttracks.s:
   hfound.Fill(seg.Plate())
#end of function, loop over files

for quarter in range (nfiles):

  condition = goodtrack+" && "+surfaces[quarter]
  trackloop(filenames[quarter],condition)

#end of loop, getting efficiency

cfound = r.TCanvas()
hfound.Draw()
cfound.Print("plots/Plate_tracks_npl15.root")
cfound.Print("plots/Plate_tracks_npl15.png")


ceff = r.TCanvas()
if (r.TEfficiency.CheckConsistency(hfound,hexpected)):
  heff = r.TEfficiency(hfound, hexpected)
  heff.Draw()
  heff.SetTitle(";Plate;#epsilon")
ceff.Print("plots/Efficiency_tracks_npl15.root")
ceff.Print("plots/Efficiency_tracks_npl15.png")

cxy = r.TCanvas()
hxy.Draw("COLZ")
cxy.Print("plots/hxy_tracks_npl15.root")
cxy.Print("plots/hxy_tracks_npl15.png")


 
