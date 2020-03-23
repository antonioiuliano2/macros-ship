''' Script to quantify how many charm daughters are tracked along the charm parent'''
import ROOT as r
import fedrarootlogon
import sys
import fedrautils
from argparse import ArgumentParser 

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

'''Used to identify extra tracks from MC'''
parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with reconstructed fedra vertices",
                    required=True)
parser.add_argument("-s", "--simulationfile", dest="simfilename", help="file with original MCtrue simulation")
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with reconstructed fedra tracks",
                    required=True)                    
parser.add_argument("-c", "--charmlist", dest="charmlistfilename", help="charm list file",
                    required=True)
parser.add_argument("-o", "--output", dest="vertexcsv", help="output csv file",
                    required=True)
options = parser.parse_args()

#IMPORTANT: PICKLE LOADS IN THE SAME ORDER OF WRITINGS. NAMES MEAN NOTHING!
with open(options.charmlistfilename, 'rb') as fp:
    charmlist = pickle.load(fp)
    daughterlist = pickle.load(fp)
    ndaughterslist = pickle.load(fp)
    decaylengthslist = pickle.load(fp)

simfile = r.TFile.Open(options.simfilename,"READ")
simtree = simfile.Get("cbmsim")

dproc = r.EdbDataProc()
gAli = dproc.PVR()
tracklist = fedrautils.buildtracks(options.tracksfilename, dproc, gAli)

#r.gSystem.Load("/home/antonio/Scrivania/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so")

vtxfile = r.TFile.Open(options.vertexfilename)
vtxtree = vtxfile.Get("vtx")

outputfile = open(options.vertexcsv,"a") 
#outputfile.write("MCEvent,charmMCTrackID,daughterMCTrackID,TrackID\n")
for itrk,track in enumerate(tracklist):
    # getting list of IDs for that MC event
    sfirst = track.GetSegmentFirst()
    MCEvent = sfirst.MCEvt()
    charmIDs = charmlist[MCEvent]
    daughterIDs = daughterlist[MCEvent]
    ndaughters = ndaughterslist[MCEvent]
    decaylengths = decaylengthslist[MCEvent]

    # getting first and last segment of the track
    slast = track.GetSegmentLast()
    #rmax = r.ShipCharmDecaySearch.FedraTrackKink(track)   
    if ((sfirst.MCTrack() in daughterIDs) and (vtxtree.GetEntries("TrackID=={}&&incoming==1".format(itrk)) == 0)):
        simtree.GetEntry(sfirst.MCEvt())
        simtracks = simtree.MCTrack #track 0 is always mother of charm
        # Topology 4: track not associated to any vertex
        outputfile.write("{0},{1},{2},{3},{4},{5},{6:.0f},{7:.0f},{8:.0f},{9},{10},{11:.0f},{12},{13:.0f},{14:.0f},{15:.0f},{16}\n".
                    format(0, 0, itrk, sfirst.MCEvt(), sfirst.MCTrack(),simtracks[0].GetPdgCode(),simtracks[0].GetPx(),simtracks[0].GetPy(),simtracks[0].GetPz(),sfirst.Aid(0),ndaughters[sfirst.Aid(0)],decaylengths[sfirst.Aid(0)] *1E+4,1, 0.,0.,0.,4))
outputfile.close()
