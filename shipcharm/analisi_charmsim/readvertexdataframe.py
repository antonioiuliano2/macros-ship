from __future__ import division
import pandas as pd
import ROOT

def df2TH1I(column, name, title, nbins,xmin,xmax ):
    '''filling a histogram from pandas'''
    histo = ROOT.TH1I(name, title, nbins, xmin, xmax)
    for item in column:
        histo.Fill(item)
    return histo

def df2TH1D(column, name, title, nbins,xmin,xmax ):
    '''filling a histogram from pandas'''
    histo = ROOT.TH1D(name, title, nbins, xmin, xmax)
    for item in column:
        histo.Fill(item)
    return histo


dfall = pd.read_csv('charmdaughtersvertices.log')

dfprimaries = dfall[dfall['topology']==1]
dfsecondaries = dfall[dfall['topology']==2]
#distance primaries vs secondaries
#dfsecondaries['MCEventID']

df = dfall[dfall['topology']!=1]


df = df.drop_duplicates(subset=["MCEventID","MCTrackID","MCMotherID"]) #many tracks are splitted up, first instance is kept

hntracks = df2TH1I(df['ntracks'],'hntracks','molteplicity of reconstructed vertices',50,0,50)
hmolt = df2TH1I(df['predmolt'],'hmolt','expected molteplicity',50,0,50)

df1 = df.groupby(['MCEventID','MCMotherID']).sum()
df1 = df1.rename(columns={"quantity":"countedmolt"})
df2 = df.groupby(['MCEventID','MCMotherID']).first() #predmolt is identical for all daughters of the same charm, I keep the first of them

#dffinal = pd.concat([df1,df2],axis=1)
dffinal = df1
dffinal['ratio'] = df1['countedmolt']/df2['predmolt']
dffinal=dffinal.fillna(0) 

entries=df1[df1['ntracks']>0].count()['countedmolt']
hratio = df2TH1D(dffinal['ratio'],'hratio','all vertices',11,0,1.1)

#cutting away primary vertices
dfcut = df[df['ntracks']<8]
df1cut = dfcut.groupby(['MCEventID','MCMotherID']).sum()
df1cut = df1cut.rename(columns={"quantity":"countedmolt"})

dfsamevertex = dfcut.groupby(['MCEventID','MCMotherID','ivtx']).count() #counting how many charm daughters are from the same FEDRA vertex
hntrackssame = df2TH1I(dfsamevertex['quantity'],'hntrackssame','n daughters from same vertex',50,0,50)

df2cut = dfcut.groupby(['MCEventID','MCMotherID']).first()

dffinalcut = df1cut

dffinalcut['ratio'] = df1cut['countedmolt']/df2cut['predmolt']
dffinalcut=dffinalcut.fillna(0) #if no daughters, we cannot reconstruct anything
entriesaftercut=df1cut[df1cut['ntracks']>0].count()['countedmolt']
hratiocut = df2TH1D(dffinalcut['ratio'],'hratiocut','no primary',11,0,1.1)

print("Larger than 0: ",entries,entriesaftercut)

#Start drawing section

cmolt = ROOT.TCanvas()
hmolt.SetLineColor(ROOT.kRed)
hmolt.SetFillColorAlpha(ROOT.kRed,0.5) # need a .rootrc file with OpenGL.CanvasPreferGL 1 to see it online
hmolt.GetXaxis().SetTitle("ndaughters")
hntracks.SetFillColorAlpha(ROOT.kBlue,0.5)
hntrackssame.SetFillColorAlpha(ROOT.kYellow,0.5)
hmolt.Draw()
hntracks.Draw("SAMES")
hntrackssame.Draw("SAMES")
cmolt.BuildLegend()

cratio = ROOT.TCanvas()
hratio.GetXaxis().SetTitle("fraction of reconstructed charm daughters")
hratio.Draw()
hratio.SetFillColorAlpha(ROOT.kBlue,0.5)
hratiocut.SetLineColor(ROOT.kRed)
hratiocut.SetFillColorAlpha(ROOT.kRed,0.5)
hratiocut.Draw("SAMES")
cratio.BuildLegend()