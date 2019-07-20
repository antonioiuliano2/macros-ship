#simple example of a script how to loop on detector points after a simulation
#syntax is 
#python -i inspectpoints.py -f inputfile.root 
#by default it looks for a file called ship.conical.Genie-TGeant4.root in the same folder
import ROOT as r
import argparse
#using an option to get the file name
def init(): #available options

  ap = argparse.ArgumentParser(
      description='Simple checks of detectors after a FairShip simulation')
  ap.add_argument('-f', '--inputfile', type=str, help="file with cbmsim simulation tree", dest='inputfile', default='ship.conical.Genie-TGeant4.root')
  args = ap.parse_args()
  return args

args = init() #getting options

#getting tree and checking if it is found
f = r.TFile(args.inputfile)
if not f:
 print 'ERROR: no file loaded. Exiting now'
 1/0
tree = f.Get('cbmsim')
if not tree:
 print 'ERROR: no cbmsim tree found. Quitting now'

nevents = tree.GetEntries() #each event is an incoming proton
print 'Number of events: ', nevents 
hmuonppt = r.TH2D("hmuonppt", "Pt vs P for muons", 400, 0, 400, 100, 0, 100)

mypdg = r.TDatabasePDG.Instance() #getting particle information
#loop on events
for ievent in range(nevents): 
 tree.GetEntry(ievent)
 tracks = tree.MCTrack #monte carlo true tracks

 for itrack, track in enumerate(tracks): #itrack provides us the TrackID in our simulation
  #getting information about tracks 
  pdgcode = track.GetPdgCode()
  procID = track.GetProcID()
  motherId = track.GetMotherId()

  startx = track.GetStartX()
  starty = track.GetStartY()
  startz = track.GetStartZ()

  px = track.GetPx()
  py = track.GetPy()
  pz = track.GetPz()
  
  momentum = track.GetP()
  transversemomentum = track.GetPt()

 # charge = 0
 # if (mypdg.GetParticle(pdgcode)): charge = mypdg.GetParticle(pdgcode).Charge()
 # else: print 'Warning, particle with pdg {} not recognized, assigned charge 0'.format(pdgcode)
  #filling histograms
  if (abs(pdgcode)==13):
   hmuonppt.Fill(momentum,transversemomentum)
#histogram setup and drawing
hmuonppt.GetXaxis().SetTitle("p[GeV]")
hmuonppt.GetYaxis().SetTitle("pt[GeV]")
hmuonppt.Draw("COLZ")
