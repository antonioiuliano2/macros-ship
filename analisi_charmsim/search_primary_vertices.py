import ROOT as r
import fedrarootlogon
from argparse import ArgumentParser 

'''Used to identify primary vertices from MC'''
parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with reconstructed fedra vertices",
                    required=True)
parser.add_argument("-o", "--output", dest="vertexcsv", help="output csv file",
                    required=True)
options = parser.parse_args()

def countX(lst, x): 
    count = 0
    for ele in lst: 
        if (ele == x): 
            count = count + 1
    return count 

def most_frequent(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = countX(List, i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num,counter 

inputfile = r.TFile.Open(options.vertexfilename) #vertexfilename
vtxtree = inputfile.Get("vtx")

dproc = r.EdbDataProc()
gAli = dproc.PVR()

outputfile = open(options.vertexcsv,'w') 

outputfile.write("ntracks,ivtx,itrk,MCEventID,MCTrackID,MCMotherID,predmolt,quantity,vx,vy,vz,topology\n") #quantity means 1 for reconstructed vertices daughter, 0 for primary or not reconstructed

for ivtx, vtx in enumerate(vtxtree):
    ntracks=vtx.n

    vID = vtx.vID
    #vertex coordinates
    vx = vtx.vx
    vy = vtx.vy
    vz = vtx.vz
    incomingtrack = most_frequent(vtx.incoming)[0] #to see if it is indeed a primary

    MCMotherIDs = vtx.MCMotherID
    MCEventIDs = vtx.MCEventID

    if (most_frequent(MCMotherIDs)[0] == -1 and incomingtrack == 1):
        #loop on tracks, to save all of them
        for itrk in range(ntracks):         
         outputfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8:.0f},{9:.0f},{10:.0f},{11}\n".format(ntracks, vID, vtx.TrackID[itrk],most_frequent(MCEventIDs)[0], vtx.MCTrackID[itrk],-1,0,0,vx,vy,vz,1))
        
outputfile.close()
# opening the file to keep only one entry for primary vertex -> the one with most tracks
import pandas as pd
df = pd.read_csv(options.vertexcsv) #outputfilename
onlyoneprimary = False
if onlyoneprimary:
 indexes = df.groupby("MCEventID")['ntracks'].idxmax()#returning indexes with maximum value
 df = df.loc[indexes] #splicing to keep only the vertices with more tracks

df.to_csv(options.vertexcsv,index=False)

