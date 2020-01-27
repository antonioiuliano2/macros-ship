import numpy as np
import pandas as pd
import sys
import ROOT

def df2TH1D(column, name, title, nbins,xmin,xmax ):
    '''filling a histogram from pandas'''
    histo = ROOT.TH1D(name, title, nbins, xmin, xmax)
    for item in column:
        histo.Fill(item)
    return histo
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

nevents = 10000
for ientry in range(0,nevents): #old style loop
    # need to know for each event how many instances of different topologies
    nbadevent = 0
    ngoodevent = 0
    nconnectedevent = 0
    nextraevent = 0
    dfprimaryvertex = dfprimaryvertices.query("MCEventID=={}".format(ientry)) #vertex for that ID
    if dfvertices.index.contains(ientry):       
        dfprimaryvertex = dfprimaryvertices.query("MCEventID=={}".format(ientry)) #vertex for that ID
        if len(dfprimaryvertex)==1:
            asprimary = dfvertices.loc[[ientry],"ivtx"]==dfprimaryvertex["ivtx"].iloc[0] #matching ivtx, they are not charm but primary vertices            
            indeces = [x[1] for x in asprimary.index]  #first index is MCEventID, second index is MCMotherID
            #need to take into account cases where only one vertex is present
            #if len(asprimary)==2:
            # indeces.append(asprimary.index[0][1])
            # indeces.append(asprimary.index[1][1])             
            #if len(asprimary)==1:
            # indeces.append(asprimary.index[0][1])
            for ndaughter in indeces:
                #print (ientry, ndaughter)
                if (asprimary.loc[ientry,ndaughter]):
                     nbad = nbad + 1
                     nbadevent = nbadevent + 1
                     #print ("Test: ",ientry,ndaughter, asprimary.iloc[ndaughter])
                else:                    
                     dfvertices.loc[(ientry,ndaughter),["topology"]] = 22 #TEST to save primary vs secondary separation
                     ngood = ngood + 1
                     ngoodevent = ngoodevent + 1
                     #print ("Test: ",ientry,ndaughter, asprimary.iloc[ndaughter])
        else: #no primary, all daughters bad
         for ndaughter in range(len(dfvertices.loc[[ientry],"ivtx"])):
             nbad = nbad +1
             nbadevent = nbadevent + 1
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
#topologymatrix = (topologymatrix + topologymatrix.T)/2.
print(topologymatrix)
print("Printing efficiencies ")
topologyefficiency = topologymatrix/nevents
print(topologyefficiency)
topologyerrors = np.sqrt(topologyefficiency*(1-topologyefficiency)/nevents) 
#(approximated) efficiency error formula

print("Printing errors")
print(topologyerrors)

#saving matrices to files
np.savetxt("MCpredictions.csv",topologymatrix,delimiter=",",fmt='%1.4e',header="Primary, Secondary, Connected, Extra, Missing")
np.savetxt("MCefficiencies.csv",topologyefficiency,delimiter=",",fmt='%1.4e',header="Primary, Secondary, Connected, Extra, Missing")
np.savetxt("MCefficiencieserrors.csv",topologyerrors,delimiter=",",fmt='%1.4e',header="Primary, Secondary, Connected, Extra, Missing")

dfprimaryvertices = dfprimaryvertices.groupby("MCEventID").first()

dftosecondary = dfvertices.query("topology==22")
dftoprimary = dfvertices.query("topology==2")

#drawing histograms
hconnectedlength = df2TH1D(dfconnected["preddecaylength"],"hconnectedlength", "connected to parent;decaylength[#mum]", 30,0.,30000. )
hextralength = df2TH1D(dfextra["preddecaylength"],"hextralength", "extra tracks;decaylength[#mum]", 30,0.,30000. )
htoprimarylength = df2TH1D(dftoprimary["preddecaylength"],"htoprimarylength", "associated to primary vertex;decaylength[#mum]", 30,0.,30000. )
htosecondarylength = df2TH1D(dftosecondary["preddecaylength"],"htosecondarylength", "associated to secondary vertex;decaylength[#mum]", 30,0.,30000. )
#normalizing
hconnectedlength.Scale(1./hconnectedlength.Integral())
hextralength.Scale(1./hextralength.Integral())
htoprimarylength.Scale(1./htoprimarylength.Integral())
htosecondarylength.Scale(1./htosecondarylength.Integral())
#setting colors
hconnectedlength.SetFillColor(ROOT.kBlue)
hextralength.SetLineColor(ROOT.kRed)
hextralength.SetFillColor(ROOT.kRed)
htoprimarylength.SetLineColor(ROOT.kYellow)
htoprimarylength.SetFillColor(ROOT.kYellow)
htosecondarylength.SetLineColor(ROOT.kMagenta)
htosecondarylength.SetFillColor(ROOT.kMagenta)
hlength = ROOT.THStack("hlength","Charm hadron decay length from true MC")
hlength.Add(htoprimarylength)
hlength.Add(htosecondarylength)
hlength.Add(hconnectedlength)
hlength.Add(hextralength)
#finally drawing
canvas = ROOT.TCanvas()
#htoprimarylength.Draw("histo")
#htosecondarylength.Draw("histo && SAMES")
#hconnectedlength.Draw("histo && SAMES")
#hextralength.Draw("histo && SAMES")
hlength.Draw("histo")
canvas.BuildLegend()

#drawing histograms
hconnectedmolt = df2TH1D(dfconnected["predmolt"],"hconnectedmolt", "connected to parent;molt[#mum]", 10,0,10 )
hextramolt = df2TH1D(dfextra["predmolt"],"hextramolt", "extra tracks;molt[#mum]", 10,0,10 )
htoprimarymolt = df2TH1D(dftoprimary["predmolt"],"htoprimarymolt", "associated to primary vertex;molt[#mum]", 10,0,10 )
htosecondarymolt = df2TH1D(dftosecondary["predmolt"],"htosecondarymolt", "associated to secondary vertex;molt[#mum]", 10,0,10 )
#normalizing
hconnectedmolt.Scale(1./hconnectedmolt.Integral())
hextramolt.Scale(1./hextramolt.Integral())
htoprimarymolt.Scale(1./htoprimarymolt.Integral())
htosecondarymolt.Scale(1./htosecondarymolt.Integral())
#setting colors
hconnectedmolt.SetFillColor(ROOT.kBlue)
hextramolt.SetLineColor(ROOT.kRed)
hextramolt.SetFillColor(ROOT.kRed)
htoprimarymolt.SetLineColor(ROOT.kYellow)
htoprimarymolt.SetFillColor(ROOT.kYellow)
htosecondarymolt.SetLineColor(ROOT.kMagenta)
htosecondarymolt.SetFillColor(ROOT.kMagenta)
hmolteplicity = ROOT.THStack("hmolteplicity","Expected molteplicity of charm decay")
hmolteplicity.Add(htoprimarymolt)
hmolteplicity.Add(htosecondarymolt)
hmolteplicity.Add(hconnectedmolt)
hmolteplicity.Add(hextramolt)
#finally drawing
canvas1 = ROOT.TCanvas()
#hextramolt.Draw("histo")
#htoprimarymolt.Draw("histo&&SAMES")
#htosecondarymolt.Draw("histo && SAMES")
#hconnectedmolt.Draw("histo && SAMES")
hmolteplicity.Draw("histo")
canvas1.BuildLegend()

#import matplotlib.pyplot as plt #using matplotlib instead of ROOT for faster drawings
#plt.figure()
#plt.hist(dfconnected["preddecaylength"],bins = 30, range = (0,30000),density=True,alpha = 0.5, label = "connected to parent")
#plt.hist(dfextra["preddecaylength"],bins = 30, range = (0,30000),density=True,alpha = 0.5,color="r", label = "extra tracks")
#plt.hist(dftoprimary["preddecaylength"],bins = 30, range = (0,30000),density=True,alpha = 0.5,color="y",label = "associated to primary vertex")
#plt.hist(dftosecondary["preddecaylength"],bins = 30, range = (0,30000),density=True,alpha = 0.5,color="m",label = "associated to secondary vertex")
#plt.xlabel("decaylength[micron]")
#plt.legend()
#plt.show()

outputlogfile = open("vertices.log","w")
selectionlogfile = open("events_withonesecondary.log","w") #list of events to check distributions on
def inspectevent(eventID,outputlogfile,selectionlogfile):
    '''inspecting the two decay topologies in the event: a vertex takes priority over connected tracks and extra tracks'''
    topologies = [-1,-1]
    outputlogfile.write("MCEvent numer {}\n".format(eventID))
    if (eventID in dfprimaryvertices.index):
            outputlogfile.write ("Found primary vertex, \n {}\n".format(dfprimaryvertices.loc[[eventID],["ivtx"]]))
    else:
        outputlogfile.write("No Primary vertex identified as such in event {}\n ".format(eventID))
    if (eventID in dftoprimary.index):        
      outputlogfile.write("Found charm daughter associated to primary,\n {}\n".format(dftoprimary.loc[[eventID],["ivtx"]]))
    if (eventID in dftosecondary.index):          
      outputlogfile.write("Charm daughter associated to secondary,\n {}\n".format(dftosecondary.loc[[eventID],["ivtx"]]))
      #look for which mother is present
      selectionlogfile.write("{}".format(eventID))
      if ((eventID,1) in dftosecondary.index):
       selectionlogfile.write(" {}".format(1))
      else:
       selectionlogfile.write(" {}".format(0)) 
      if ((eventID,2) in dftosecondary.index):
       selectionlogfile.write(" {}".format(2))
      else:
       selectionlogfile.write(" {}".format(0))          
     # for single in dfvertices.loc[eventID]:
       # help(single)
      selectionlogfile.write("\n")
    if (eventID in dfconnected.index):
      outputlogfile.write("Track Connected to parent, \n {}\n".format(dfconnected.loc[[eventID],["itrk"]]))
    if (eventID in dfextra.index):
      outputlogfile.write("Extra track, \n {}\n".format(dfextra.loc[[eventID],["itrk"]]))

#dfvertices.groupby(['MCEventID','ivtx']).sum()
# Double Signal: topology 1, 2 ,2 total 5
# Single Signal: topology 1,2 or 2,2: total less than 5

print("Saving vertex log info")
for ievent in range(nevents):
  inspectevent(ievent,outputlogfile,selectionlogfile)
outputlogfile.close()
selectionlogfile.close()
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
