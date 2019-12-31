import ROOT as r
import fedrarootlogon
from argparse import ArgumentParser 

'''Used to identify secondary vertices from MC'''
parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with reconstructed fedra vertices",
                    required=True)
parser.add_argument("-c", "--charmlist", dest="charmlistfilename", help="charm list file",
                    required=True)
parser.add_argument("-o", "--output", dest="vertexcsv", help="output csv file",
                    required=True)
options = parser.parse_args()
# python search_secondary_vertices.py inputvertexfile outputcsv inputcharm
try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

newversion = True

#IMPORTANT: PICKLE LOADS IN THE SAME ORDER OF WRITINGS. NAMES MEAN NOTHING! (see writecharmdaughters.py to check)
with open(options.charmlistfilename, 'rb') as fp:
    charmlist = pickle.load(fp)
    daughterlist = pickle.load(fp)
    ndaughterslist = pickle.load(fp)
    decaylengthlist = pickle.load(fp)

def countX(lst, x): 
    count = 0
    for ele in lst: 
        if (ele == x): 
            count = count + 1
    return count 

def most_frequent(List): 
    '''find the most frequent instance in a list'''
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = countX(List, i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num,counter

inputfile = r.TFile.Open(options.vertexfilename,"READ")
vtxtree = inputfile.Get("vtx")

outputfile = open(options.vertexcsv,"a") 
#main loop on vertices
#outputfile.write("ntracks,ivtx,itrk,MCEventID,MCTrackID,MCMotherID,predmolt,preddecaylen,quantity,vx,vy,vz,topology\n") #quantity means 1 for reconstructed vertices daughter, 0 for not reconstructed
for ivtx,vtx in enumerate(vtxtree):
    ntracks=vtx.n

    vID = vtx.vID
    vx = vtx.vx
    vy = vtx.vy
    vz = vtx.vz

    segments = vtx.s
    incominglist = vtx.incoming # 1 if track starts from there
    #loop on vertices associated to track
    nprevioussegments = 0
    for (MCEventID,MCTrackID,MCMotherID,TrackID,nseg,incomingtrack) in zip(vtx.MCEventID,vtx.MCTrackID,vtx.MCMotherID,vtx.TrackID,vtx.nseg,incominglist):
        #first and last segment
        firstsegment = segments[0+nprevioussegments]
        lastsegment = segments[nseg-1+nprevioussegments]
        nprevioussegments = nprevioussegments + nseg
                
        charmIDs = charmlist[MCEventID]
        daughterIDs = daughterlist[MCEventID]
        ndaughters = ndaughterslist[MCEventID]
        decaylengths = decaylengthlist[MCEventID]  # in cm, needs to be multiplied for 1E+4 to micron conversion   
        if ((firstsegment.MCTrack() in charmIDs) and (lastsegment.MCTrack() in daughterIDs)):
        # Topology 3: track connected to parent
         outputfile.write("{0},{1},{2},{3},{4},{5},{6},{7:.0f},{8},{9:.0f},{10:.0f},{11:.0f},{12}\n".\
                    format(ntracks, vID, TrackID, lastsegment.MCEvt(), lastsegment.MCTrack(),lastsegment.Aid(0),ndaughters[lastsegment.Aid(0)],decaylengths[lastsegment.Aid(0)] *1E+4,1, vx,vy,vz,3))
        #getting list of IDs for charm daughters from that event  
        #if MCMotherID == -1 and incomingtrack == 1:
        #        outputfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8:.0f},{9:.0f},{10:.0f},{11}\n".format(ntracks, vID, TrackID, MCEventID, MCTrackID,MCMotherID,0,1,vx,vy,vz,1))
        if MCTrackID in daughterIDs and incomingtrack == 1:
                outputfile.write("{0},{1},{2},{3},{4},{5},{6},{7:.0f},{8},{9:.0f},{10:.0f},{11:.0f},{12}\n".\
                    format(ntracks, vID, TrackID, MCEventID, MCTrackID,MCMotherID,ndaughters[MCMotherID],decaylengths[MCMotherID] *1E+4,1, vx,vy,vz,2))

def lookforcharm(df,ievent,charmID):
    ''' handling the case where no charm daughter was reconstructed, missing entry in dataframe'''
    checkdf = df.query("MCEventID == "+str(ievent)+" and MCMotherID == "+str(charmID))
    if checkdf.empty: return False
    else: return True
def addmissing(inputfile):
    import pandas as pd
    df = pd.read_csv(inputfile,header=0) #header = 1 if there is a dummy line due to Load Fedra libs, sep is needed if it is not a true csv file
    #we need to add zero for missing charm daughters
    outputfile = open(inputfile,"a")  

    nevents = len(charmlist)
    for ievent in range(nevents):
        charmIDs = charmlist[ievent]
        ndaughters = ndaughterslist[ievent]
        for ID in charmIDs: #looping on the two charmIDs
            found = lookforcharm(df,ievent,ID)        
            if not found: outputfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8:.0f},{9:.0f},{10:.0f},{11}\n".format(0, 0, 0, ievent, 0, ID, ndaughters[ID], 0,0.,0.,0.,0))

    outputfile.close()

outputfile.close()
#addmissing('charmdaughtersvertices.log')
