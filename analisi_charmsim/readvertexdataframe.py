from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('charmdaughtersvertices.xlsx')

df = df.drop_duplicates(subset=["MCEventID","MCTrackID"]) #many tracks are splitted up, first instance is kept
plt.figure()
plt.hist(df['ntracks'],50,[0,50],alpha=0.5,label='molteplicity of reconstructed vertices')
plt.legend()
df1 = df.groupby(['MCEventID','MCMotherID']).count()
#remove now useless columns
df1.pop('ntracks')
df1.pop('ivtx')
df1.pop('itrk')
df1.pop('MCTrackID')
df1 = df1.rename(columns={"molt":"countedmolt"})

df2 = df.drop_duplicates(subset=['MCEventID','MCMotherID']) #to compare with predicted
df2 = df.groupby(['MCEventID','MCMotherID']).max()

df2.pop('ntracks')
df2.pop('ivtx')
df2.pop('itrk')
df2.pop('MCTrackID')

dffinal = pd.concat([df1,df2],axis=1)

dffinal['ratio'] = dffinal['countedmolt']/dffinal['molt']

entries=df1.count()['countedmolt']
plt.figure()
plt.hist(dffinal['ratio'],10,[0,1],label='all vertices: '+str(entries)+' charm decays',alpha=0.5)

#cutting away primary vertices
dfcut = df[df['ntracks']<=df['molt']]
df1cut = dfcut.groupby(['MCEventID','MCMotherID']).count()
#remove now useless columns
df1cut.pop('ntracks')
df1cut.pop('ivtx')
df1cut.pop('itrk')
df1cut.pop('MCTrackID')
df1cut = df1cut.rename(columns={"molt":"countedmolt"})

df2cut = dfcut.drop_duplicates(subset=['MCEventID','MCMotherID']) #to compare with predicted
df2cut = dfcut.groupby(['MCEventID','MCMotherID']).max()

df2cut.pop('ntracks')
df2cut.pop('ivtx')
df2cut.pop('itrk')
df2cut.pop('MCTrackID')

dffinalcut = pd.concat([df1cut,df2cut],axis=1)

dffinalcut['ratio'] = dffinalcut['countedmolt']/dffinalcut['molt']
entriesaftercut=df1cut.count()['countedmolt']
plt.hist(dffinalcut['ratio'],10,[0,1],label='no primary '+str(entriesaftercut)+' charm decays',alpha=0.5)

print(entries,entriesaftercut)
plt.legend()
plt.xlabel('Fraction reconstructed daughters')
plt.show()