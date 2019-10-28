import ROOT as r
import recognizecharmdaughters
from argparse import ArgumentParser 

'''Used to save a dictionary to associate event with charm daugters'''
parser = ArgumentParser()
parser.add_argument("-s", "--inputsim", dest="simulationfilename", help="file with cbmsim simulation",
                    required=True)
parser.add_argument("-co", "--charmlistoutput", dest="charmlistfilename", help="output charm list file",
                    required=True)
options = parser.parse_args()

inputfile = r.TFile.Open(options.simulationfilename)
inputtree = inputfile.Get("cbmsim")

nentries = inputtree.GetEntries()
charmlist = []
daughterlist = []
ndaughterslist = []
decaylengthlist = []

for ievent in range(nentries):
 charmIDs, daughterIDs, ndaughters, decaylengths = recognizecharmdaughters.getdaughtertracks(inputtree,ievent)

 charmlist.append(charmIDs)
 daughterlist.append(daughterIDs)
 ndaughterslist.append(ndaughters)
 decaylengthlist.append(decaylengths)

print ("Saving lists")

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

with open(options.charmlistfilename, 'wb') as fp:
#IMPORTANT: PICKLE LOADS IN THE SAME ORDER OF WRITINGS. NAMES MEAN NOTHING!
    pickle.dump(charmlist, fp)
    pickle.dump(daughterlist, fp)
    pickle.dump(ndaughterslist,fp)
    pickle.dump(decaylengthlist,fp)