import ROOT as r
import fedrarootlogon
import sys

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

inputfile = r.TFile.Open(sys.argv[1],"READ")
vertexrec = inputfile.Get("EdbVertexRec")
#IMPORTANT: PICKLE LOADS IN THE SAME ORDER OF WRITINGS. NAMES MEAN NOTHING!
with open('charmlist.p', 'rb') as fp:
    charmlist = pickle.load(fp)
    daughterlist = pickle.load(fp)
    ndaughterslist = pickle.load(fp)

vtxlist = vertexrec.eVTX
#main loop on vertices
print ("ntracks", "ivtx", "itrk", "MCEventID", "MCTrackID","MCMotherID")
for ivtx,vtx in enumerate(vtxlist):
    ntracks=vtx.N()
    #loop on vertices associated to track
    for itrk in range(ntracks):
        track = vtx.GetTrack(itrk)
        motherid = track.Aid(0)
        MCTrackID = track.MCTrack()
        MCEventID = track.MCEvt()
        #getting list of IDs for charm daughters from that event
        daughterIDs = daughterlist[MCEventID]
        ndaughters = ndaughterslist[MCEventID]
        if MCTrackID in daughterIDs:
                print (ntracks, ivtx, itrk, MCEventID, MCTrackID,motherid,ndaughters[motherid])

