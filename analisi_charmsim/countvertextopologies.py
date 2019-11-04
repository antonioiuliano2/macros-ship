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

dfconnected = dfcharm[dfcharm['topology'] > 2.5]

print ("How many connected? ",len(dfconnected))

dfextra = dfcharm[dfcharm['topology'] > 3.5]

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
for ientry in range(1000): #old style loop
    if dfcharm.index.contains(ientry):
        dfprimaryvertex = dfprimaryvertices.query("MCEventID=={}".format(ientry)) #vertex for that ID
        if len(dfprimaryvertex)==1:
            dfcharm.loc[[ientry],"ivtx"]==dfprimaryvertex["ivtx"].iloc[0] #matching ivtx, they are not charm but primary vertices


dfvertices = dfcharm[dfcharm['topology']<=2.5]
dfvertices = dfvertices[dfvertices['topology']>=1.5]

dfprimaryvertices["MCEventID"]==dfvertices["MCEventID"]
#dfvertices.groupby(['MCEventID','ivtx']).sum()
# Double Signal: topology 1, 2 ,2 total 5
# Single Signal: topology 1,2 or 2,2: total less than 5

#Quick check, molteplicity (TO DO; uso topology)
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
