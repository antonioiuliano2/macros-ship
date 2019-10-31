import pandas as pd

#getting the dataframe

dfall = pd.read_csv('MC_vertexlist_testremainingtracks.csv')

#removing duplicates of same track

dfall = dfall.sort_values(["MCEventID","topology"])

df = dfall.drop_duplicates(subset=["MCEventID","MCTrackID","MCMotherID"]) #many tracks are splitted up, first instance is kept

dfcharm = df.groupby(['MCEventID','MCMotherID']).first()

dfconnected = dfcharm[dfcharm['topology'] > 2.5]

print ("How many connected? ",len(dfconnected))

dfextra = dfcharm[dfcharm['topology'] > 3.5]

print ("How many extra? ",len(dfextra))

pairconnected = dfconnected.groupby(['MCEventID']).sum()
print (len(pairconnected))
print (len(pairconnected[pairconnected['quantity']==2]))

#Vertices
dfvertices = dfcharm[dfcharm['topology'] < 2.5]
dfvertices = dfvertices[dfvertices['topology']>1.5]
print("How many vertices? ", len(dfvertices))
print(len(dfvertices))
dfprimary = dfvertices[dfvertices['ntracks']>=6]
print("Associated to primary? ", len(dfprimary))
dfsecondary = dfvertices[dfvertices['ntracks']<6]
print("Reconstructed as secondary? ", len(dfsecondary))