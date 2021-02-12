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
nmotheridtrue = truesimdf.groupby(["MCEvent","MotherID"]).count()

simdf_vertices = simdf.query("VertexS >= 0")

#taking only one entry per track

simdf_vertices =  simdf_vertices.sort_values(["MCEvent", "MCTrack", "PID"],ascending = [True,True,False]) #from upstream to downstream plates
simdf_vertices =  simdf_vertices.groupby(["MCEvent","MCTrack"]).first().reset_index()

print("size of track sample with vertices", len(simdf))

#sorting per vertexID

simdf_vertices = simdf_vertices.sort_values("VertexS")

#how many tracks per vertex
ntracks = simdf_vertices.groupby("VertexS").count()["VertexE"].to_numpy()

#what is the most frequent MCmotherId per vertex, and how many tracks have it? (I am assuming same motherID from differents events are not present)
mostfrequentmotherId = simdf_vertices.groupby("VertexS")["MotherID"].apply(lambda x: x.value_counts().index[0])
ntracks_mostfrequentmotherId = simdf_vertices.groupby("VertexS")["MotherID"].apply(lambda x: x.value_counts().values[0])

#only one entry per vertex
vertexdf = simdf_vertices.groupby("VertexS").first()

#adding found branches
vertexdf["ntracks"] = ntracks
vertexdf["ntracks_true"] = ntracks_mostfrequentmotherId
vertexdf["VertexMCID"] = mostfrequentmotherId

#purity

vertexdf["purity"] = vertexdf["ntracks_true"] / vertexdf["ntracks"]