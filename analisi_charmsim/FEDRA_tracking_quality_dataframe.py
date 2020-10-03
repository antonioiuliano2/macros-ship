'''evaluating tracking goodness and purity with respect to the MC true track'''
from __future__ import division
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import seaborn as sns


simdf = pd.read_csv(sys.argv[1])

print("original size of segment sample ", len(simdf))
nsegtrue = simdf.groupby(["MCEvent","MCTrack"]).count()


#keeping only segments which have been connected to tracks, sorting the dataframe
simdf = simdf.query("TrackID>=0")
print("size of tracked segment sample ", len(simdf))

simdf = simdf.sort_values(["MCEvent", "MCTrack", "PID"],ascending = [True,True,False]) #from upstream to downstream plates

def efficiency(simdf,nsegtrue):
 '''first, group by MC particle, this tell us how many segments belong to the same particle'''

 #then, require them to belong to the same reconstructed track

 nsegreco = simdf.groupby(["MCEvent","MCTrack","TrackID"]).count()

 #keeping only the first segment, then adding the number of segments as new column
 simdftracks = simdf.groupby(["MCEvent","MCTrack","TrackID"]).first()
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
 nseg = simdf.groupby("TrackID").count()
 recotrackdf = simdf.groupby("TrackID").first()
 recotrackdf["nseg"] = nseg["ID"]

 nsegreco = simdf.groupby(["MCEvent","MCTrack","TrackID"]).count()
 #keeping only the first segment, then adding the number of segments as new column
 simdftracks = simdf.groupby(["MCEvent","MCTrack","TrackID"]).first()
 simdftracks["nsegreco"] = nsegreco["ID"]

 #sorting with descending nsegreco, so we keep the track with largest nsegreco
 simdftracks.sort_values(["MCEvent","MCTrack","nsegreco"],ascending=[True,True,False])
 simdftracks = simdftracks.reset_index()
 simdfbesttrack = simdftracks.groupby("TrackID").first()
 simdfbesttrack["nseg"] = nseg["ID"]
 simdfbesttrack["purity"] = simdfbesttrack["nsegreco"] / simdfbesttrack["nseg"]

 return simdfbesttrack

simdfbesttrack = efficiency(simdf,nsegtrue)
simdfbesttrack2 = purity(simdf)

#profile histogram, computing mean of the column for each segments
length = np.linspace(1,29,29)
efficiencies = []
for l in length:
    subset = simdfbesttrack.query("nsegtrue == {}".format(l))
    efficiencies.append(subset.mean()["efficiency"])

plt.plot(length,efficiencies,"r")
plt.show()