import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()
# <code to time>

X_68_all = np.zeros(0)
Y_68_all = np.zeros(0)
W_68_all = np.zeros(0)

xrootd_prefix = "root://eospublic.cern.ch/"

path = "/eos/experiment/ship/user/edursov/pycondor_out/muons_in_ms_sc_unrolled_spill_full_info_for_warm_part/11/"

#loop over files
for ientry in range(0,10):
 print(ientry)
 simfile = uproot.open(xrootd_prefix+path+"{}/ship.conical.MuonBack-TGeant4.root".format(ientry))
 cbmsim = simfile["cbmsim"]

 #positions
 X_68 = cbmsim["sco_68Point"]["sco_68Point.fX"].array()
 Y_68 = cbmsim["sco_68Point"]["sco_68Point.fY"].array()
 #Z_68 = cbmsim["sco_68Point"]["sco_68Point.fZ"].array()
 #PdgCode_68 = cbmsim["sco_68Point"]["sco_68Point.fPdgCode"].array()
 TrackID_68 = cbmsim["sco_68Point"]["sco_68Point.fTrackID"].array()
 W = cbmsim["MCTrack"]["MCTrack.fW"].array()

 #need to add weight from MCTrack to arrays, with TrackID info
 #just find where is not empty (eventID) and the track value for that hit (TrackID)
 nhits=ak.num(TrackID_68)
 eventIDs = np.where(nhits>0)
 TrackIDs = TrackID_68[nhits>0]

 W_68 = W[eventIDs][TrackIDs]

 #concatenate array for each entry
 X_68_all = np.concatenate([X_68_all,ak.flatten(X_68).to_numpy()])
 Y_68_all = np.concatenate([Y_68_all,ak.flatten(Y_68).to_numpy()])
 W_68_all = np.concatenate([W_68_all,ak.flatten(W_68).to_numpy()])

end = time.time()

print(f"Time taken to run the code was {end-start} seconds")
#plot x and y
plt.figure()
plt.hist2d(X_68_all,Y_68_all,weights = W_68_all,range=[[-100,100],[-100,100]],bins=[200,200],cmin=1)
plt.xlabel("x[cm]")
plt.ylabel("y[cm]")
plt.colorbar()
plt.show()