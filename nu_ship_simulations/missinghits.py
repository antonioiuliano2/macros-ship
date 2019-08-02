'''Are there missing TargetPoints for a given track?(29 July 2019)'''
import ROOT as r
import ctypes #to pass as reference to C classes
import sys

inputfile = r.TFile.Open(sys.argv[1])
inputtree = inputfile.Get("cbmsim")

nentries = inputtree.GetEntries()

target = r.Target()

def dumphits(nevent,ntrack):
 inputtree.GetEntry(nevent)
 targetpoints = inputtree.TargetPoint
 for point in targetpoints:

     NWall = ctypes.c_int(0)
     NRow = ctypes.c_int(0)
     NColumn = ctypes.c_int(0)
     NPlate = ctypes.c_int(0)
     EmCES = ctypes.c_bool(False)
     EmBrick = ctypes.c_bool(False)
     EmTop = ctypes.c_bool(False)

     if (point.GetTrackID() == ntrack):
         detID = point.GetDetectorID()
         target.DecodeBrickID(detID,NWall,NRow,NColumn,NPlate,EmCES,EmBrick,EmTop)
         print (detID, NWall, NPlate, point.GetZ(), point.GetEnergyLoss())
print "(detID, NWall, NPlate, Zposition, Eloss)"
dumphits(4,1)