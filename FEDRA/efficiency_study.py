#python version of efficiency_studyC, to easily hande multiple files
import ROOT as r
import fedrarootlogon

goodtrack ="t.eFlag>=0 &&t.eProb>0.01&&npl >= 15"

surfaces = {}

surfaces[0] = "t.eX < 62500 && t.eY<50000"
surfaces[1] = "t.eX > 62500 && t.eY<50000"
surfaces[2] = "t.eX < 62500 && t.eY>50000"
surfaces[3] = "t.eX > 62500 && t.eY>50000"

filenames = {}
filenames[0] = "linked_tracks_firstquarter.root"
filenames[1] = "linked_tracks_secondquarter.root"
filenames[2] = "linked_tracks_thirdquarter.root"
filenames[3] = "linked_tracks_fourthquarter.root"

#degining histograms
hexpected = r.TH1D("hexpected", "Tracks expected to be found in each plate", 58,0,58)
hfound = r.TH1D("hfound", "Tracks with an associated segment in each plate", 58,0,58); 

hxy = r.TH2D("hxy", "2D position of tracks;x[mm];y[mm]",120,0,120,100,0,100);

nquarters = 4

for quarter in range (nquarters):
 #getting file and tree

 firstfile = r.TFile.Open(filenames[quarter])
 firsttracks = firstfile.Get("tracks")

 #reading entries from selection
 firsttracks.Draw(">>lst",r.TCut(goodtrack+" && "+surfaces[quarter]))

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
  hexpected.Fill(0)
  #next, found plates
  for seg in firsttracks.s:
   hfound.Fill(seg.Plate())
#end of loop, getting efficiency
if (r.TEfficiency.CheckConsistency(hfound,hexpected)):
  heff = r.TEfficiency(hfound, hexpected);
  heff.Draw();
  heff.SetTitle(";npl;#epsilon");

cxy = r.TCanvas()
hxy.Draw("COLZ")

 
