import ROOT
import sys
'''Trees from SciFi Event Builder have different entry numbers. But they can be identified according to the main timestamp. How much time sensitivity? '''
globalfile = ROOT.TFile.Open(sys.argv[1])
scififile = ROOT.TFile.Open(sys.argv[2])

#outputfile = ROOT.TFile.Open("entrycomparison.root","UPDATE")

globaltree = globalfile.Get("cbmsim")
scifitree = scififile.Get("cbmsim")

sensitivity = 1e+6 #10 ms
#start loop on first tree
hentries = ROOT.TH1I((ROOT.TString("hentries_")+ROOT.TString(sys.argv[3])).Data(),"How many entries for event time in {:.2e} s".format(sensitivity),5,0,5)

for event in globaltree:
 globaltime = globaltree.EventHeader.GetEventTime()
 #we need to check if only one entry if the other tree matches this condition
 nentries = scifitree.GetEntries("TMath::Abs(EventHeader.GetEventTime()-"+str(globaltime)+")<"+str(sensitivity))
 hentries.Fill(nentries)

hdelta = ROOT.TH1D("hdelta", "histogram of time differences", 400,-2,2)
for entry in range(scifitree.GetEntries()):
    scifitree.GetEntry(entry)
    globaltree.GetEntry(entry)
    globalheader = scifitree.EventHeader
    scifiheader = scifitree.EventHeader
    globalheader = globaltree.EventHeader
    delta = globalheader.GetEventTime() - scifiheader.GetEventTime()
    hdelta.Fill(delta/1E+6)

hdelta.GetXaxis().SetTitle("#Deltat[ms]")
hdelta.SetTitle("globalheadereventtime - scifiheadereventtime")
ROOT.gStyle.SetOptStat(111111)
hdelta.Draw()
#outputfile.Write()
