'''script to draw volumes from ROOT geometry
   usage: python -i drawrootvolume.py geofile.root
   Example of volumes names: tTauNuDet, Wall, volTarget, MuonShieldArea, volMuFilter, etc.
'''
import ROOT as r
import sys
#usage: python -i drawrootvolume.py geofile.root
r.TGeoManager.Import(sys.argv[1])
#r.gGeoManager.GetTopVolume().Draw("ogl")
#r.gGeoManager.GetVolume("MuonShieldArea").Draw("ogl")
r.gGeoManager.GetVolume("MuonDetector").Draw("ogl")
r.gGeoManager.SetVisLevel(9)
