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

total_entries = 0
total_W = 0

hxy = r.TH2D("hxy","check muon rates;x[m];y[m]", 40,-0.2,0.2,40,-0.2,0.2)

for f in os.listdir(files_dir):
    if not f.endswith('pkl'): continue
    print("File name", f)
    print("Numer entries", get_data(files_dir+"/"+f).shape[1])
    px,py,pz,x,y,z,particle,W = get_data(files_dir+"/"+f)
    
    for x_i, y_i, W_i in zip( x, y, W):
     hxy.Fill(x_i, y_i, W_i)
    
    total_entries = total_entries + get_data(files_dir+"/"+f).shape[1]
    total_W = total_W + np.sum(W)

#plotting the xy data

cxy = r.TCanvas()
hxy.Draw("COLZ")

print("Total entries", total_entries)
print("Total weight", total_W)

print("Weighted muons in acceptance", hxy.Integral())
