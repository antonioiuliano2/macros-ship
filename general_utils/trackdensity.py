'''just a quick code to evaluate track density within all configurations'''

import ROOT as r
import fedrarootlogon
import sys

cxy = r.TCanvas()
cpid = r.TCanvas()


nplates=57
surfaces = []

#need to use splicing to remove default first argument
filenames = sys.argv[1:]
nfiles = len(filenames)
print(nfiles)
if nfiles==1:
 surfaces.append("1")
elif nfiles==4:
 surfaces.append("t.eX > 0.     && t.eX < 62500.  && t.eY > 0.   && t.eY<50000. ")
 surfaces.append("t.eX > 62500. && t.eX < 125000. && t.eY > 0.   && t.eY<50000. ")
 surfaces.append("t.eX > 0.     && t.eX < 62500.  && t.eY>50000.  && t.eY<100000.")
 surfaces.append("t.eX > 62500. && t.eX < 125000. && t.eY>50000.  && t.eY<100000.")

 #surfaces.append("t.eX > 0.     && t.eX < 62500.  && t.eY > 41000.   && t.eY<50000. ")
 #surfaces.append("t.eX > 62500. && t.eX < 125000. && t.eY > 41000.   && t.eY<50000. ")
 #surfaces.append("t.eX > 0.     && t.eX < 62500.  && t.eY>50000.  && t.eY<51000.")
 #surfaces.append("t.eX > 62500. && t.eX < 125000. && t.eY>50000.  && t.eY<51000.")

#histogram for track xy distribution
hxytot = r.TH2D("hxytot","XY distribution of tracks on the surface;x[#mum];y[#mum]",125,0,125000,100,0,100000)
hpidtot = r.TH1D("hpidtot","Plate of track segments;PID;tracks/cm^2",nplates,0,nplates)

iquarter = 0
for filename,surface in zip(filenames,surfaces):

 trackfile = r.TFile.Open(filename,"READ")
 tracktree = trackfile.Get("tracks")

 #draw xy
 tracktree.Draw("t.eY:t.eX>>hxy{}(125,0,125000,100,0,100000)".format(iquarter),surface)
 hxy = r.gDirectory.Get("hxy{}".format(iquarter))

 #draw pid
 tracktree.Draw("s.ePID>>hpid{}({},0,{})".format(iquarter,nplates,nplates),surface)
 hpid = r.gDirectory.Get("hpid{}".format(iquarter))

 hxytot.Add(hxy)
 iquarter = iquarter+1
 hpidtot.Add(hpid)
#drawing found histograms
cxy.cd()
hxytot.Draw("COLZ")
cpid.cd()
hpidtot.Scale(1./125)
hpidtot.Draw("histo")
