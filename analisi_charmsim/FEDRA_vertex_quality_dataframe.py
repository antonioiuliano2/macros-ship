'''evaluating vertex reconstruction goodness and purity with respect to the MC mother ID'''
from __future__ import division
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

simdf = pd.read_csv(sys.argv[1])

print("original size of segment sample ", len(simdf))

#counting how many tracks are with same event and mother (all dataset, before selection)
truesimdf = simdf.sort_values(["MCEvent","MCTrack"])
truesimdf = truesimdf.groupby(["MCEvent","MCTrack"]).first()
nmotheridtrue = truesimdf.groupby(["MCEvent","MotherID"]).count()["ID"].to_numpy()

#datasetwithonlyothers
trackmothers = truesimdf.groupby(["MCEvent","MotherID"]).first()
trackmothers["ndaughters"] = nmotheridtrue
trackmothers = trackmothers.reset_index()

simdf_vertices = simdf.query("VertexS >= 0")

#taking only one entry per track

simdf_vertices =  simdf_vertices.sort_values(["MCEvent", "FEDRATrackID", "PID"],ascending = [True,True,False]) #from upstream to downstream plates
simdf_vertices =  simdf_vertices.groupby(["MCEvent","FEDRATrackID"]).first().reset_index()

print("size of track sample with vertices", len(truesimdf))

#sorting per vertexID

simdf_vertices = simdf_vertices.sort_values("VertexS")

#how many tracks per vertex
ntracks = simdf_vertices.groupby("VertexS").count()["VertexE"].to_numpy()

print("grouping vertices and estimating purity")
#what is the most frequent MCmotherId per vertex, and how many tracks have it? (I am assuming same motherID from differents events are not present)
mostfrequentmotherId = simdf_vertices.groupby("VertexS")["MotherID"].apply(lambda x: x.value_counts().index[0])
ntracks_mostfrequentmotherId = simdf_vertices.groupby("VertexS")["MotherID"].apply(lambda x: x.value_counts().values[0])

mostfrequentProcessId = simdf_vertices.groupby("VertexS")["ProcID"].apply(lambda x: x.value_counts().index[0])

#only one entry per vertex
vertexdf = simdf_vertices.groupby("VertexS").first()

#adding found branches
vertexdf["ntracks"] = ntracks
vertexdf["ntracks_true"] = ntracks_mostfrequentmotherId
vertexdf["VertexMCID"] = mostfrequentmotherId
#renaming motherid to match the new definition
del vertexdf["MotherID"]
vertexdf["MotherID"] = vertexdf["VertexMCID"]
del vertexdf["ProcID"]
vertexdf["ProcID"] = mostfrequentProcessId

#only the most abundant vertex for MCID
vertexdf.sort_values(["VertexMCID","ntracks"],ascending=[True,False])
vertexdf.groupby("VertexMCID").first()
#purity

vertexdf["purity"] = vertexdf["ntracks_true"] / vertexdf["ntracks"]

mergeddf = vertexdf.merge(trackmothers, how = 'left', on=["MCEvent","MotherID"])

mergeddf["efficiency"] = mergeddf["ntracks_true"]/mergeddf["ndaughters"]

plt.figure()
plt.hist(mergeddf["purity"],range = [0,1], color = 'r',label = 'purity')
plt.hist(mergeddf["efficiency"],range=[0,1], label = 'efficiency')
plt.xlabel("Percentage (%)",fontsize=20)
plt.legend(fontsize=20)

#only hadronic interactions
plt.figure()
hadintdf = mergeddf.query("ProcID_x == 23")
plt.hist(hadintdf["purity"],range = [0,1], color = 'r')
plt.hist(hadintdf["efficiency"],range=[0,1])
plt.xlabel("Percentage (%)",fontsize=20)
plt.legend(fontsize=20)
plt.show()



