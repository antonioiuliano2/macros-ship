import uproot
import sys
import matplotlib.pyplot as plt
import numpy as np

shipcharmtree = uproot.open(sys.argv[1])["shippositions"]

timestamp = shipcharmtree.array("eventtime")

plt.plot(timestamp,"b*",label = "data from a spill of RUN 2814")
plt.xlabel("ientry",fontsize="large")
plt.ylabel("timestamp[ns]",fontsize="large")

nentries = len(timestamp)
counter = np.linspace(1,nentries,nentries)

coefficients = np.polyfit(counter,timestamp,1)

estimatedtimestamp = counter * coefficients[0] + coefficients[1]

print ("fit results: ",coefficients[0], coefficients[1])

plt.plot(estimatedtimestamp,"r",label="{:.2E}".format(coefficients[0])+"*x+"+"{:.2E}".format(coefficients[1]))
plt.legend()
plt.show()

#import ROOT 

#graph = ROOT.TGraph()
#for i, time in enumerate(timestamp):
  #  graph.SetPoint(i,i,time)
#graph.Draw("AP*")
#ROOT.gStyle.SetOptFit(1111)
#graph.Fit("pol1")