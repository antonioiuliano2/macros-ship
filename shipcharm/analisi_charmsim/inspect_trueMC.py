'''user script to inspect properties of tracks and event'''
import ROOT as r
import sys

inputfile = r.TFile.Open(sys.argv[1])
simtree = inputfile.Get("cbmsim")

pdgdatabase = r.TDatabasePDG.Instance()

def inspecttrack(MCEventID, MCTrackID):
    '''printing information about track within event'''
    simtree.GetEntry(MCEventID)
    tracks = simtree.MCTrack

    selectedtrack = tracks[MCTrackID]

    print("Track PDG: ",selectedtrack.GetPdgCode()," momentum ",selectedtrack.GetP(), " produced by ", selectedtrack.GetProcName())
    if (selectedtrack.GetMotherId() > -1 ):
     print("Mother PDG: ",tracks[selectedtrack.GetMotherId()].GetPdgCode()," momentum ",tracks[selectedtrack.GetMotherId()].GetP())

def inspectlist(MCEventID, mytracks):

    '''printing information about list of tracks within an event'''
    simtree.GetEntry(MCEventID)  
    tracks = simtree.MCTrack

    for MCTrackID in mytracks:
      selectedtrack = tracks[MCTrackID]

      print("Track PDG: ",selectedtrack.GetPdgCode()," momentum ",selectedtrack.GetP(), " produced by ", selectedtrack.GetProcName())
      if (selectedtrack.GetMotherId() > -1 ):
       print("Mother PDG: ",tracks[selectedtrack.GetMotherId()].GetPdgCode()," momentum ",tracks[selectedtrack.GetMotherId()].GetP())
      print("\n")
