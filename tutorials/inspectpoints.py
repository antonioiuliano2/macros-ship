#simple example of a script how to loop on detector points after a simulation
#syntax is 
#python -i inspectpoints.py -f inputfile.root

import ROOT as r
import argparse
#using an option to get the file name
def init(): #available options

  ap = argparse.ArgumentParser(
      description='Simple checks of detectors after a FairShip simulation')
  ap.add_argument('-f', '--inputfile', type=str, help="file with cbmsim simulation tree", dest='inputfile', default='None')
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
hbeamxy = r.TH2D("hbeamxy", "Simulated profile of the beam", 200, -1 ,1, 200, -1, 1)
hxy = r.TH2D('hxy', '2D distribution for first Rpc plane',200,-100,100,150,-50,100)

for ievent in range(nevents): #loop on events
 tree.GetEntry(ievent)
 tracks = tree.MCTrack #monte carlo true tracks
 #getting information about tracks 

 startx = tracks[0].GetStartX()
 starty = tracks[0].GetStartY()
 hbeamxy.Fill(startx, starty)
 rpcpoint = tree.MuonTaggerPoint
 for hit in rpcpoint: #loop on MuonTaggerPoint
  pdgcode = hit.PdgCode()
  if (r.TMath.Abs(pdgcode) == 2112) or (pdgcode == 22):
   continue #not physical neutral hits  

  #getting information about hits
  detID = hit.GetDetectorID()  
  trackID = hit.GetTrackID()
  nplane = (detID - 4000)/10000

  x = hit.GetX()
  y = hit.GetY()
  z = hit.GetZ()
  px = hit.GetPx()
  py = hit.GetPy()
  pz = hit.GetPz()
  
  momentum = pow(pow(px, 2) + pow(py, 2) + pow(pz,2),0.5)
  #filling histograms
  if (nplane == 1):
   hxy.Fill(x,y)

outputfile = r.TFile("quickcheck.root","RECREATE") #save canvas to a file for later checks
c1 = r.TCanvas()
hbeamxy.Draw()
hbeamxy.GetXaxis().SetName('x[cm]')
hbeamxy.GetYaxis().SetName('y[cm]')
c2 = r.TCanvas()
hxy.Draw() 
hxy.GetXaxis().SetName('x[cm]')
hxy.GetYaxis().SetName('y[cm]')

hbeamxy.Write()
hxy.Write()
