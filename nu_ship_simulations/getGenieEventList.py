#I want a list of Genie Events in target acceptance

import ROOT as r

#cbmsim->Scan("GenieEventID:MCEventHeader.fEventId:MCTrack[0].GetEnergy()")

inputfile = r.TFile.Open("inECC_ship.conical.Genie-TGeant4.root")
inputsim = inputfile.Get("cbmsim")

outfile = open("ineccgenieevtnumber_numubarccdis.txt","w") #to read genie events in acceptance

#start event loop
for entry in inputsim:

 header = entry.MCEventHeader
 ttreenumber = entry.GenieEventID #to be updated
 origmceventid = header.GetEventID() 

 genieeventid = (ttreenumber-2) * 1000 + (origmceventid-1)

 outfile.write("{}\n".format(genieeventid))

outfile.close()



