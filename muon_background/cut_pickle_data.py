import gzip
import pickle
import numpy as np
import ROOT as r
import os
import matplotlib.pyplot as plt

import pandas as pd

files_dir = '/eos/user/l/lcattela/muons/files_middle' #Directory with files
def get_data(file):
    with gzip.open(file,'rb') as f:
        data = pickle.load(f)
    return data

PDG = r.TDatabasePDG.Instance()
m_muon = PDG.GetParticle(13).Mass()
m_proton = PDG.GetParticle(2212).Mass()

#all_data = []
nbins = 120

total_entries = 0
total_W = 0

df_acc = pd.DataFrame()

totalpcounts = np.zeros(400,dtype=int)
totalcounts = np.zeros(nbins,dtype = int) #we cannot concatenate all data due to memory, so we bin them for each file
for f in os.listdir(files_dir):
    if not f.endswith('pkl'): continue
    print("File name", f)
    print("Numer entries", get_data(files_dir+"/"+f).shape[1])
    px,py,pz,x,y,z,particle,W = get_data(files_dir+"/"+f)
    
    df = pd.DataFrame({"px":px,"py":py,"pz":pz,"x":x,"y":y,"z":z,"particle":particle,"W":W})
    
    df_cut = df.query("abs(x)<0.2 and abs(y)<0.2")
    
    df_acc = pd.concat([df_acc,df_cut])
    
#all_data = np.concatenate(all_data,axis = 1)
#writing muons in acceptance into output file
outputfile = gzip.open("muons_acceptance_middle.pkl","wb")
pickle.dump([df_acc["px"].to_numpy(),df_acc["py"].to_numpy(),df_acc["pz"].to_numpy(),
             df_acc["x"].to_numpy(),df_acc["y"].to_numpy(),df_acc["z"].to_numpy(),
             df_acc["particle"].to_numpy(),df_acc["W"].to_numpy()],outputfile)
outputfile.close()

#px,py,pz,x,y,z,particle,W = all_data #momentum, position [in meters], particle_id [13 or -13], muons weight
#Pt = np.sqrt(px**2+py**2)
#P = np.sqrt(Pt**2+pz**2)
#plotting Ecm

plt.figure()
plt.hist2d(df_acc["x"],df_acc["y"],bins=[100,100],range=[[-0.2,0.2],[-0.2,0.2]])
plt.xlabel("x[m]")
plt.ylabel("y[m]")
plt.show()

