#Allows to load FEDRA libraries in PyROOT
# Before runnning this script do:
# cp src/*/*.pcm   lib
# cp src/*/*/*.pcm lib

from ROOT import gSystem
gSystem.Load('libMatrix.so');
gSystem.Load('libTree.so');
gSystem.Load('libHist.so');
gSystem.Load('libPhysics.so');
gSystem.Load('libEve.so');
gSystem.Load('libGeom.so');
gSystem.Load('libEve.so');

gSystem.Load('libvt.so');
gSystem.Load('libEphys.so');
gSystem.Load('libEmath.so');
gSystem.Load('libEdb.so');
gSystem.Load('libEbase.so');
gSystem.Load('libEdr.so');
gSystem.Load('libEIO.so');
gSystem.Load('libAlignment.so');
gSystem.Load('libAnalysis.so');
gSystem.Load('libScan.so');
gSystem.Load('libDataConversion.so');

gSystem.Load('libEGA.so');
gSystem.Load('libEdd.so');
gSystem.Load('libEMC.so');
gSystem.Load('libShower.so');
gSystem.Load('libEmr.so');
gSystem.Load('libEDA.so');

print ("Load FEDRA libs")
