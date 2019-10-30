import pandas as pd

#getting the dataframe

dfall = pd.read_csv('MC_vertexlist_testremainingtracks.csv')

#removing duplicates of same track

df = dfall.drop_duplicates(subset=["MCEventID","MCTrackID","MCMotherID"]) #many tracks are splitted up, first instance is kept

dfcharm = df.groupby(['MCEventID','MCMotherID']).first()

dfconnected = dfcharm[dfcharm['topology'] > 2.5]

print ("How many connected? ",len(dfconnected))

dfextra = dfcharm[dfcharm['topology'] > 3.5]

print ("How many extra? ",len(dfextra))

testpairconnected = dfconnected.groupby(['MCEventID']).sum()
print (len(testpairconnected))
print (len(testpairconnected[testpairconnected['quantity']==2]))