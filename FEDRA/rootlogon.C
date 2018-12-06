// Before runnning this script do:
// cp src/*/*.pcm   lib
// cp src/*/*/*.pcm lib

R__LOAD_LIBRARY(libMatrix.so);
R__LOAD_LIBRARY(libTree.so);
R__LOAD_LIBRARY(libHist.so);
R__LOAD_LIBRARY(libPhysics.so);
R__LOAD_LIBRARY(libEve.so);
R__LOAD_LIBRARY(libGeom.so);
R__LOAD_LIBRARY(libEve.so);

R__LOAD_LIBRARY(libvt.so);
R__LOAD_LIBRARY(libEphys.so);
R__LOAD_LIBRARY(libEmath.so);
R__LOAD_LIBRARY(libEdb.so);
R__LOAD_LIBRARY(libEbase.so);
R__LOAD_LIBRARY(libEdr.so);
R__LOAD_LIBRARY(libEIO.so);
R__LOAD_LIBRARY(libAlignment.so);
R__LOAD_LIBRARY(libAnalysis.so);
R__LOAD_LIBRARY(libScan.so);
R__LOAD_LIBRARY(libDataConversion.so);

R__LOAD_LIBRARY(libEGA.so);
R__LOAD_LIBRARY(libEdd.so);
R__LOAD_LIBRARY(libEMC.so);
R__LOAD_LIBRARY(libShower.so);
R__LOAD_LIBRARY(libEmr.so);
R__LOAD_LIBRARY(libEDA.so);

void rootlogon()
{
    cout << "Load FEDRA libs" << endl;
}