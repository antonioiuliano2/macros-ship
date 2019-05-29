#transcription of some methods used by Giuliana for her decay search
import ROOT
import fedrarootlogon

def Kinkangle(parenttracktx, parenttrackty, daughtertrackstxlist, daughtertrackstylist):
 '''computing the average kink angle, between mother and daughter particles'''
 kink=0.

 parenttx = parenttracktx
 parentty = parenttrackty 

 for daughtertx, daughterty in zip(daughtertrackstxlist, daughtertrackstylist):

  daughtertx = daughtertrack.TX()
  daughterty = daughtertrack.TY()
  kink += sqrt(pow(parenttx - daughtertx,2) + pow(parentty - daughterty,2))
 
 ndaughters = len(daughtertrackslist)
 return kink/ndaughters

def GetpT(parenttracktx, parenttrackty, daughtertracktx, daughtertrackty, daughtermomentum):
 '''transverse momentum, of daughter with respect with mother particle'''
 parenttx = parenttracktx
 parentty = parenttrackty
 daughtertx = daughtertracktx
 daughterty = daughtertrackty
 
 costh = (1+parenttx*daughtertx + parentty*daughterty)/sqrt((1+parenttx*parenttx+parentty*parentty)*(1+daughtertx*daughtertx+daughterty*daughterty))

 pT2 = daughtermomentum * ROOT.TMath.Sin(ROOT.TMath.ACos(costh))
 
 return pT2


def IPtoVertex(vertexpos, trackstartpos, tracktx, trackty):
 '''derived from Elena IPtoTarget method, three TVector3 as inputs. Remember the b'''
 #getting vertex and start of track position as two TVector3 objects
 #vertexpos = ROOT.TVector3(vx,vy,vz)

 print "Start of Track: ",trackstartpos(0), trackstartpos(1), trackstartpos(2)
 print "Vertex position: ",vertexpos(0), vertexpos(1), vertexpos(2)
 
 dz = vertexpos(2) - trackstartpos(2)
 ipx = tracktx * dz + trackstartpos(0)
 ipy = trackty * dz + trackstartpos(1)

 ip = ROOT.TMath.Sqrt(pow(ipx,2)+pow(ipy,2))
