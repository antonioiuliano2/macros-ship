'''print geometry information according to loaded fairship geometry'''
import ROOT as r
from ShipGeoConfig import ConfigRegistry
globalDesigns = {'2016':{'dy':10.,'dv':5,'ds':7,'nud':1,'caloDesign':0,'strawDesign':4},\
                 '2018':{'dy':10.,'dv':6,'ds':9,'nud':3,'caloDesign':3,'strawDesign':10}}
default = '2018'

dy           = globalDesigns[default]['dy'] # max height of vacuum tank
dv           = globalDesigns[default]['dv'] # 4=TP elliptical tank design, 5 = optimized conical rectangular design, 6=5 without segment-1
ds           = globalDesigns[default]['ds'] # 5=TP muon shield, 6=magnetized hadron, 7=short magnet design, 9=optimised with T4 as constraint, 8=requires config file
                                            # 10=with field map for hadron absorber
nud          = globalDesigns[default]['nud'] # 0=TP, 1=new magnet option for short muon shield, 2= no magnet surrounding neutrino detector
caloDesign   = globalDesigns[default]['caloDesign'] # 0=ECAL/HCAL TP  1=ECAL/HCAL TP + preshower 2=splitCal  3=ECAL/ passive HCAL 
strawDesign  = globalDesigns[default]['strawDesign'] # simplistic tracker design,  4=sophisticated straw tube design, horizontal wires (default), 10=2cm straw diameter for 2018 layout
geofile = None

ship_geo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py", Yheight = dy, tankDesign = dv, \
                                                muShieldDesign = ds, nuTauTargetDesign=nud, CaloDesign=caloDesign, strawDesign=strawDesign, muShieldGeo=geofile)

#ship_geo.NuTauTarget
print("Start loop on EmuMagnet parameters")
for attr, value in ship_geo.EmuMagnet.__dict__.items():
        print (attr, value)
print("Start loop on NuTauTarget parameters")
for attr, value in ship_geo.NuTauTarget.__dict__.items():
        print (attr, value)
print("Start loop on NuTauTT parameters")
for attr, value in ship_geo.NuTauTT.__dict__.items():
        print (attr, value)
print("Start loop on tauHPT parameters")
for attr, value in ship_geo.tauHPT.__dict__.items():
        print (attr,value)
print("Start loop on tauMudet parameters")
for attr, value in ship_geo.tauMudet.__dict__.items():
        print (attr, value)
