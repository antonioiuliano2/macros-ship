import numpy as np
import pandas as pd
import sys
#getting the dataframe

dfall = pd.read_csv(sys.argv[1])

dfall = dfall.sort_values(["MCEventID","topology","ntracks"])
#removing duplicates of same track

dfall = dfall.sort_values(["MCEventID","topology"])

df = dfall.drop_duplicates(subset=["MCEventID","MCTrackID","MCMotherID"]) #many tracks are splitted up, first instance is kept

dfcharm = df[df['topology']>1.5].groupby(['MCEventID','MCMotherID']).first()

print ("How many charm daughters in log?", len(dfcharm))

dfconnected = dfcharm[dfcharm['topology'] == 3]

print ("How many connected? ",len(dfconnected))

dfextra = dfcharm[dfcharm['topology'] == 4]

print ("How many extra? ",len(dfextra))

#counting for pair
pairconnected = dfconnected.groupby(['MCEventID']).sum()
print ("Single charm daughter connected to parent: ",len(pairconnected))
print ("Both charm daughters connected to parent: ",len(pairconnected[pairconnected['quantity']==2]))

pairextra = dfextra.groupby(['MCEventID']).sum()
print ("Single charm daughter connected to parent: ",len(pairextra))
print ("Both charm daughters connected to parent: ",len(pairextra[pairextra['quantity']==2]))

#vertex reconstruction (1 primary, 2 secondary)

#recognizing primary vertices
dfprimaryvertices =df[df["topology"]==1]
indexes = dfprimaryvertices.groupby("MCEventID")['ntracks'].idxmax()#returning indexes with maximum value
dfprimaryvertices = dfprimaryvertices.loc[indexes] #splicing to keep only the vertices with more tracks
#come cerco i vertici primari di un dato evento da quelli secondari?
nbad = 0
ngood = 0

dfvertices = dfcharm[dfcharm['topology']==2]

topologymatrix = np.zeros([5,5])

for ientry in range(1000): #old style loop
    # need to know for each event how many instances of different topologies
    nbadevent = 0
    ngoodevent = 0
    nconnectedevent = 0
    nextraevent = 0
    if dfvertices.index.contains(ientry):       
        dfprimaryvertex = dfprimaryvertices.query("MCEventID=={}".format(ientry)) #vertex for that ID
        if len(dfprimaryvertex)==1:
            asprimary = dfvertices.loc[[ientry],"ivtx"]==dfprimaryvertex["ivtx"].iloc[0] #matching ivtx, they are not charm but primary vertices
            for ndaughter in range(len(asprimary)):
                if (asprimary.iloc[ndaughter]):
                     nbad = nbad + 1
                     nbadevent = nbadevent + 1
                     #print ("Test: ",ientry,ndaughter, asprimary.iloc[ndaughter])
                else:
                     ngood = ngood + 1
                     ngoodevent = ngoodevent + 1
                     #print ("Test: ",ientry,ndaughter, asprimary.iloc[ndaughter])
        else: #no primary, all daughters bad
         for ndaughter in range(len(dfvertices.loc[[ientry],"ivtx"])):
             nbad = nbad +1
    # how many connected tracks and extra tracks?
    if (ientry in dfconnected.index):
      nconnectedevent = len(dfconnected.loc[ientry])
    if (ientry in dfextra.index):
      nextraevent = len(dfextra.loc[ientry])
    tottopologiesevent = nbadevent + ngoodevent +nextraevent+ nconnectedevent
    if (tottopologiesevent > 2):
      print ("ERROR: TOO MANY TOPOLOGIES!")
      1./0.
    nmissingevent = 2 - tottopologiesevent
    if (ientry %100 ==0): 
      print("Test for event {}: {} as primary, {} as secondary, {} connected to parent, {} extratracks, {} missingevent".format(ientry,nbadevent,ngoodevent,nconnectedevent,nextraevent,nmissingevent))
    #matrix filling: starting from diagonal
    if nbadevent == 2:
     topologymatrix[0,0] = topologymatrix[0,0] + 1       
    if ngoodevent == 2:
     topologymatrix[1,1] = topologymatrix[1,1] + 1
    if nconnectedevent == 2:
     topologymatrix[2,2] = topologymatrix[2,2] + 1
    if nextraevent == 2:
     topologymatrix[3,3] = topologymatrix[3,3] + 1
    if nmissingevent == 2:
     topologymatrix[4,4] = topologymatrix[4,4] + 1   
    #off diagonal, first row:
    if nbadevent == 1 and ngoodevent == 1:
      topologymatrix[0,1] = topologymatrix[0,1] + 1   
    if nbadevent == 1 and nconnectedevent == 1:
      topologymatrix[0,2] = topologymatrix[0,2] + 1   
    if nbadevent == 1 and nextraevent == 1:
      topologymatrix[0,3] = topologymatrix[0,3] + 1       
    if nbadevent == 1 and nmissingevent == 1:
      topologymatrix[0,4] = topologymatrix[0,4] + 1   
    #off diagonal, second row
    if ngoodevent == 1 and nconnectedevent == 1:
     topologymatrix[1,2] = topologymatrix[1,2] + 1
    if ngoodevent == 1 and nextraevent == 1:
     topologymatrix[1,3] = topologymatrix[1,3] + 1
    if ngoodevent == 1 and nmissingevent == 1:
     topologymatrix[1,4] = topologymatrix[1,4] + 1
    #off diagonal, third row
    if nconnectedevent == 1 and nextraevent == 1:
     topologymatrix[2,3] = topologymatrix[2,3] + 1
    if nconnectedevent == 1 and nmissingevent == 1:
     topologymatrix[2,4] = topologymatrix[2,4] + 1
    #off diagonal, fourth row
    if nextraevent == 1 and nmissingevent == 1:
     topologymatrix[3,4] = topologymatrix[3,4] + 1    
#end of loop
print("After loop: found {} as secondary, {} as primary".format(ngood,nbad))
print(np.sum(topologymatrix))
#obtaining a symmetric matrix
topologymatrix = (topologymatrix + topologymatrix.T)/2.
print(topologymatrix)
dfprimaryvertices = dfprimaryvertices.groupby("MCEventID").first()


def inspectevent(eventID):
    '''inspecting the two decay topologies in the event: a vertex takes priority over connected tracks and extra tracks'''
    topologies = [-1,-1]
    if (eventID in dfprimaryvertices.index):
            print ("Primary vertex")
            print (dfprimaryvertices.loc[[eventID],["ivtx"]])
    else:
        print ("No Primary vertex identified as such in event ",eventID)
    if (eventID in dfvertices.index):    
      nsecondaries = len(dfvertices.loc[eventID])
      print ("How many secondary vertices", nsecondaries)
      #print (dfvertices.loc[[eventID],["ivtx"]])
     # for single in dfvertices.loc[eventID]:
       # help(single)
    if (eventID in dfconnected.index):
      print ("Track Connected to parent")
      print (dfconnected.loc[[eventID],["itrk"]])
    if (eventID in dfextra.index):
      print ("Extra track")
      print (dfextra.loc[[eventID],["itrk"]])

#dfvertices.groupby(['MCEventID','ivtx']).sum()
# Double Signal: topology 1, 2 ,2 total 5
# Single Signal: topology 1,2 or 2,2: total less than 5


def rawcheck():
 '''Quick check, molteplicity (now replaced with topology)'''
 dfprimary = dfvertices[dfvertices['ntracks']>=6]
 print ("How many associated to primary",len(dfprimary))
 dfsecondary = dfvertices[dfvertices['ntracks']<6]
 print ("How many associated to secondary",len(dfsecondary))

 #counting for pair
 pairprimary = dfprimary.groupby(['MCEventID']).sum()
 print ("Single charm daughter connected to primary: ",len(pairprimary))
 print ("Both charm daughters connected to primary: ",len(pairprimary[pairprimary['quantity']==2])) 

 pairsecondary = dfsecondary.groupby(['MCEventID']).sum()
 print ("Single charm daughter connected to secondary: ",len(pairsecondary))
 print ("Both charm daughters connected to secondary: ",len(pairsecondary[pairsecondary['quantity']==2]))
