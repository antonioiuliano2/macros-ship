#remake tracks, root reading issue by Valerio

import ROOT as r
import fedrarootlogon
import fedrautils
from argparse import ArgumentParser 

parser = ArgumentParser()
parser.add_argument("-o", "--output", dest="outputfilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)


options = parser.parse_args()
#reading input tracks
dproc = r.EdbDataProc()
gAli = dproc.PVR()
inputtracklist = fedrautils.buildtracks(options.tracksfilename, dproc, gAli)

#writing output tracks
xv = 0.
yv = 0.
dproc.MakeTracksTree(inputtracklist, xv,yv,options.outputfilename)
