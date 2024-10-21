import gzip
import pickle
import numpy as np
import os
files_dir = 'root:://eosuser.cern.ch/eos/user/l/lcattela/muons/files' #Directory with files
def get_data(file):
    with gzip.open(file,'rb') as f:
        data = pickle.load(f)
    return data

all_data = []
for f in os.listdir(files_dir):
    if not f.endswith('pkl'): continue
    all_data.append(get_data(f))
all_data = np.concatenate(all_data,axis = 1)

px,py,pz,x,y,z,particle,W = all_data #momentum, position [in meters], particle_id [13 or -13], muons weight
Pt = np.sqrt(px**2+py**2)
P = np.sqrt(Pt**2+pz**2)
