'''evaluating tracking goodness and purity with respect to the MC true track'''
from __future__ import division
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


simdf = pd.read_csv(sys.argv[1])

if ("FEDRATrackID" not in simdf.keys()):
 simdf["FEDRATrackID"] = simdf["TrackID"]

if ("quarter" not in simdf.keys()):
 simdf["quarter"] = 0

simdf["P"] = np.sqrt(np.power(simdf["startPx"],2) + np.power(simdf["startPy"],2) + np.power(simdf["startPz"],2)) #trimomentum
print("original size of segment sample ", len(simdf))
nsegtrue = simdf.groupby(["MCEvent","MCTrack"]).count()


#keeping only segments which have been connected to tracks, sorting the dataframe
simdf = simdf.query("FEDRATrackID>=0")
print("size of tracked segment sample ", len(simdf))

simdf = simdf.sort_values(["MCEvent", "MCTrack", "PID"],ascending = [True,True,False]) #from upstream to downstream plates

def efficiency(simdf,nsegtrue):
 '''first, group by MC particle, this tell us how many segments belong to the same particle'''

 #then, require them to belong to the same reconstructed track

 nsegreco = simdf.groupby(["MCEvent","MCTrack","FEDRATrackID","quarter"]).count()

 #keeping only the first segment, then adding the number of segments as new column
 simdftracks = simdf.groupby(["MCEvent","MCTrack","FEDRATrackID","quarter"]).first()
 simdftracks["nsegreco"] = nsegreco["ID"]

 #sorting with descending nsegreco, so we keep the track with largest nsegreco
 simdftracks.sort_values(["MCEvent","MCTrack","nsegreco"],ascending=[True,True,False])
 simdftracks = simdftracks.reset_index()
 simdfbesttrack = simdftracks.groupby(["MCEvent","MCTrack"]).first()

 simdfbesttrack["nsegtrue"] = nsegtrue["ID"]
 simdfbesttrack["efficiency"] = simdfbesttrack["nsegreco"] / simdfbesttrack["nsegtrue"]

 return simdfbesttrack

def purity(simdf):
 '''for each reco track, divide number of segments belonging to the same MCTrack for the total number of reconstructed segments'''
 #first, count nseg for each reconstructed track
 nseg = simdf.groupby(["FEDRATrackID","quarter"]).count()
 recotrackdf = simdf.groupby(["FEDRATrackID","quarter"]).first()
 recotrackdf["nseg"] = nseg["ID"]

 nsegreco = simdf.groupby(["MCEvent","MCTrack","FEDRATrackID","quarter"]).count()
 #keeping only the first segment, then adding the number of segments as new column
 simdftracks = simdf.groupby(["MCEvent","MCTrack","FEDRATrackID","quarter"]).first()
 simdftracks["nsegreco"] = nsegreco["ID"]

 #sorting with descending nsegreco, so we keep the track with largest nsegreco
 simdftracks.sort_values(["MCEvent","MCTrack","nsegreco"],ascending=[True,True,False])
 simdftracks = simdftracks.reset_index()
 simdfbesttrack = simdftracks.groupby(["FEDRATrackID","quarter"]).first()
 simdfbesttrack["nseg"] = nseg["ID"]
 simdfbesttrack["purity"] = simdfbesttrack["nsegreco"] / simdfbesttrack["nseg"]
 
 return simdfbesttrack

#obtaining and merging samples

simdfbesttrack = efficiency(simdf,nsegtrue)
print("Computed efficiency")
simdfbesttrack2 = purity(simdf)
print("Computed purity")

#mergedsimdf = simdfbesttrack2.merge(simdfbesttrack)

#profile histogram, computing mean of the column for each segments
minlength = 3 #useless count 1 and 2 segments
length = np.linspace(minlength,29,30-minlength)
efficiencies = []
purities = []

efficiencies_errors = []
purities_errors = []
for l in length:

    print("Adding tracks with {} segments".format(l))
    subset = simdfbesttrack.query("nsegtrue == {}".format(l))
    subset2 = simdfbesttrack2.query("nsegreco == {}".format(l))
    N = len(subset)
    N2 = len(subset2)
    
    efficiencies.append(subset.mean()["efficiency"]*100)
    efficiencies_errors.append(subset.std()["efficiency"]/np.sqrt(N)*100)

    purities.append(subset2.mean()["purity"]*100)
    purities_errors.append(subset2.std()["purity"]/np.sqrt(N2)*100)

#making efficiency/purity vs P

efficiencies_momentum = []
purities_momentum = []

efficiencies_errors_momentum = []
purities_errors_momentum = []

minP = 0
maxP = 400
Pstep = 20
Pbins = np.arange(minP, maxP, Pstep)

for P in Pbins:
    print("Adding tracks with P between {} and {}".format(P,P+Pstep))
    
    subset = simdfbesttrack.query("P >= {} and P < {}".format(P, P + Pstep))
    subset2 = simdfbesttrack2.query("P >= {} and P < {}".format(P, P + Pstep))

    N = len(subset)
    N2 = len(subset2)
    
    efficiencies_momentum.append(subset.mean()["efficiency"]*100)
    efficiencies_errors_momentum.append(subset.std()["efficiency"]/np.sqrt(N)*100)

    purities_momentum.append(subset2.mean()["purity"]*100)
    purities_errors_momentum.append(subset2.std()["purity"]/np.sqrt(N2)*100)

#drawing figures
plt.figure()
plt.errorbar(length,efficiencies,c = "r",fmt='.', yerr=efficiencies_errors,ecolor="r",label="efficiency",markersize=18)
plt.errorbar(length,purities,c = "b",fmt='.',yerr=purities_errors,ecolor="b",label= "purity",markersize=18)
plt.xlabel("Track length (number of segments)",fontsize=20)
plt.ylabel("Percentage (%)",fontsize=20)
plt.ylim((0, 100)) 
plt.legend(fontsize=20)

plt.figure()
plt.errorbar(Pbins,efficiencies_momentum,c = "r",fmt='.', yerr=efficiencies_errors_momentum,ecolor="r",label="efficiency",markersize=18)
plt.errorbar(Pbins,purities_momentum,c = "b",fmt='.',yerr=purities_errors_momentum,ecolor="b",label= "purity",markersize=18)
plt.xlabel("Momentum [GeV/c]",fontsize=20)
plt.ylabel("Percentage (%)",fontsize=20)
plt.ylim((0, 100)) 
plt.legend(fontsize=20)
plt.show()
