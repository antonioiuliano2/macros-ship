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

maxndaughters = 100
#branches to write in a tree
charmIDs = array( 'i', 2 * [ 0 ] )
ndaughters1 = array ('i', 1 * [0] )
ndaughters2 = array ('i', 1* [0])
daughter1IDs = array( 'i', maxndaughters*[ 0. ] )
daughter2IDs = array( 'i', maxndaughters*[0. ] )

outputfile = TFile( 'charmids.root', 'recreate' )
charminfo = TTree( 'charminfo', 'tree with charm and daughters MCIDs' )

#setting branches
charminfo.Branch( 'charmIDs', charmIDs, 'charmIDs[2]/I')
charminfo.Branch( 'ndaughters1', ndaughters1, 'ndaughters1/I' )
charminfo.Branch( 'ndaughters2', ndaughters2, 'ndaughters2/I' )
charminfo.Branch( 'daughter1IDs', daughter1IDs, 'daughter1IDs[ndaughters1]/I' )
charminfo.Branch( 'daughter2IDs', daughter2IDs, 'daughter2IDs[ndaughters2]/I' )


for ievent in range(nentries):
 charmIDs, daughterIDs, ndaughters, decaylengths = recognizecharmdaughters.getdaughtertracks(inputtree,ievent)

 charmlist.append(charmIDs)
 daughterlist.append(daughterIDs)
 ndaughterslist.append(ndaughters)
 decaylengthlist.append(decaylengths)

 #values to pass to tree
 #charm ids
 charmIDs[0] = charmlist[0]
 charmIDs[1] = charmlist[1]
 #ndaughters
 ndaughters1 = ndaughters[0]
 ndaughters2 = ndaughters[1]
 #daughtersIDs
 for idaughter in range(ndaughters1):
  daughters1IDs[idaughter] = daughterIDs[idaughter]

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
