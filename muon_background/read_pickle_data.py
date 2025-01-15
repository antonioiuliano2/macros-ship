import gzip
import pickle
import numpy as np
import ROOT as r
import os
import matplotlib.pyplot as plt
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


totalpcounts = np.zeros(400,dtype=int)
totalcounts = np.zeros(nbins,dtype = int) #we cannot concatenate all data due to memory, so we bin them for each file
for f in os.listdir(files_dir):
    if not f.endswith('pkl'): continue
    print("File name", f)
    print("Numer entries", get_data(files_dir+"/"+f).shape[1])
    px,py,pz,x,y,z,particle,W = get_data(files_dir+"/"+f)
    
    total_entries = total_entries + get_data(files_dir+"/"+f).shape[1]
    total_W = total_W + np.sum(W)
    #compute ECM (28 GeV is the 400GeV proton over proton maximum for SHiP, as reference)
    p = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    E = np.sqrt(m_muon**2 + p ** 2)
    Ecm = np.sqrt(2*E*m_proton)
    #print("Min ECM ",np.min(Ecm),"Min P", np.min(p))

    counts, bins = np.histogram(Ecm,bins=nbins,range=[0,30])
    pcounts, bins = np.histogram(p,bins=400,range=[0,400])
    
    totalcounts = totalcounts + counts
    totalpcounts = totalpcounts + pcounts
    #all_data.append(get_data(f))
#all_data = np.concatenate(all_data,axis = 1)

#px,py,pz,x,y,z,particle,W = all_data #momentum, position [in meters], particle_id [13 or -13], muons weight
#Pt = np.sqrt(px**2+py**2)
#P = np.sqrt(Pt**2+pz**2)
#plotting Ecm

print("Total entries", total_entries)
print("Total weight", total_W)

