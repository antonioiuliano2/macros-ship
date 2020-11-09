import ROOT,os,sys
import __builtin__ as builtin
import rootUtils as ut
import shipunit as u
import shipRoot_conf

shipRoot_conf.configure()
#getting the file and test points
inputfile = ROOT.TFile.Open("ship.conical.Genie-TGeant4.root")
inputtree = inputfile.Get("cbmsim")

fgeo = ROOT.TFile.Open("geofile_full.conical.Genie-TGeant4.root")

#load geometry
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')

import shipDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for creating VMC field
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)
# run.Init()
gMan = fgeo.FAIRGeom
import geomGeant4

if hasattr(ShipGeo.Bfield,"fieldMap"):
  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True,withVirtualMC = False)

print("DONE PREPARATION")

#init geometry and mag.field

#gMan = ROOT.gGeoManager #probably already instantiated in fgeo.FAIRGeom



bfield = ROOT.genfit.FairShipFields()
print("GETTING GLOBAL FIELD")

bfield.setField(fieldMaker.getGlobalField())

fM = ROOT.genfit.FieldManager.getInstance()
print("FIELD INIT")
fM.init(bfield)

geoMat = ROOT.genfit.TGeoMaterialInterface()

ROOT.genfit.MaterialEffects.getInstance().init(geoMat)

#init fitter
  #fitter          = ROOT.genfit.KalmanFitter()
  #fitter          = ROOT.genfit.KalmanFitterRefTrack()
fitter = ROOT.genfit.DAF()
fitter.setMaxIterations(50)


inputtree.GetEntry(1)
testpoints = []
#adding muon' hits manually
testpoints.append(inputtree.HptPoint[18])
testpoints.append(inputtree.HptPoint[20])
testpoints.append(inputtree.HptPoint[22])
testpoints.append(inputtree.HptPoint[24])
testpoints.append(inputtree.HptPoint[26])

def findTracks():
    pdg    = testpoints[0].PdgCode()
 #   if not self.PDG.GetParticle(pdg): continue # unknown particle   

    posM = ROOT.TVector3(0, 0, 0)
    momM = ROOT.TVector3(0,0,3.*u.GeV)
# approximate covariance
    covM = ROOT.TMatrixDSym(6)
    resolution = 0.

 # trackrep
    rep = ROOT.genfit.RKTrackRep(pdg)
# smeared start state
    stateSmeared = ROOT.genfit.MeasuredStateOnPlane(rep)
    rep.setPosMomCov(stateSmeared, posM, momM, covM)
# create track
    seedState = ROOT.TVectorD(6)
    seedCov   = ROOT.TMatrixDSym(6)
    rep.get6DStateCov(stateSmeared, seedState, seedCov)
    theTrack = ROOT.genfit.Track(rep, seedState, seedCov)
    hitCov = ROOT.TMatrixDSym(3)
    #hitCov[6][6] = resolution*resolution
    for point in testpoints:
      tp = ROOT.genfit.TrackPoint(theTrack) # note how the point is told which track it belongs to 

      hitpos = ROOT.TVectorD(3)
      hitpos[0] = point.GetX()
      hitpos[1] = point.GetY()
      hitpos[2] = point.GetZ()
  #SpacepointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
      measurement = ROOT.genfit.SpacepointMeasurement(hitpos, hitCov, 1, 6, tp)   

#   measurement = ROOT.genfit.WireMeasurement(m,hitCov,1,6,tp) # the measurement is told which trackpoint it belongs to
      # print measurement.getMaxDistance()
      measurement.setMaxDistance(ShipGeo.strawtubes.InnerStrawDiameter/2.)
      # measurement.setLeftRightResolution(-1)
      tp.addRawMeasurement(measurement) # package measurement in the TrackPoint                                          
      theTrack.insertPoint(tp)  # add point to Track
   # print "debug meas",atrack,nM,stationCrossed[atrack],self.sTree.MCTrack[atrack],pdg
   #trackCandidates.append([theTrack,atrack])
    if not theTrack.checkConsistency():
     print('Problem with track before fit, not consistent')
    else:
     fitter.processTrack(theTrack)
     fitStatus   = theTrack.getFitStatus()
     nmeas = fitStatus.getNdf()   
     chi2        = fitStatus.getChi2()/nmeas 
     print("Chi square: ", chi2)

findTracks()
