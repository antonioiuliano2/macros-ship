'''first, we need to find track ID'''
import fedrarootlogon
import ROOT as r
import sys

tracksfile = r.TFile.Open("/home/antonio/Dottorato/CharmData/decay_search_MC/CH1_charmcascade_23_12_19/secondquarter/linked_tracks.root")

trackstree = tracksfile.Get("tracks")

trackID = -1
trackstree.BuildIndex("s[0].ePID","s[0].eID")
#trackstree.SetBranchAddress("trid",trackID)

def readtree(xstart, ystart):
 '''find track IDs generating showers in cell starting at xstart and ystart (size 1 cm2)'''
 print("Arrived at file",xstart,ystart)
 showerfile = r.TFile.Open("/home/antonio/cernbox/Synched/Charmsim_Showreco/Showreco_ds_25_01/Shower_{}_{}.root".format(xstart,ystart))
 showertree = showerfile.Get("treebranch")

 outputfile = open("/home/antonio/cernbox/Synched/Charmsim_Showreco/Showreco_ds_25_01/showerlog_{}_{}.dat".format(xstart,ystart),"w")
 #loop on found showers
 outputfile.write("IShower TrackID StarterID StarterPlate Sizeb Output15 Output30\n")
 for ishower, shower in enumerate(showertree):
    segmentsid = showertree.idb
    segmentsplate = showertree.plateb
    segmentsz = showertree.zb
    #I am interested in the selector segment, which was the first segment of a track

    firstsegz = min(segmentsz)
    for starterz, starterid,starterplate in zip(segmentsz, segmentsid, segmentsplate):
     if (starterz == firstsegz): #find selector segment(s)

      trackID = -1
     #getting the track information
      pippo = trackstree.GetEntryWithIndex(starterplate,starterid)

      trackID = trackstree.GetEntryNumberWithIndex(starterplate,starterid)

      #NN output
      if (trackID > -1):
       outputfile.write("{} {} {} {} {} {} {}\n".format(ishower, trackID, starterid,starterplate, showertree.sizeb, showertree.output15,showertree.output30))

 outputfile.close()

for xcode in range(6,12):
    for ycode in range(5):
        xstart = xcode * 10000
        ystart = ycode * 10000
        readtree(xstart, ystart)
        

    