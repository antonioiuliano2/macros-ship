#test tracking and reconstruction with MC information

import ROOT as r
import fedrarootlogon
from collections import Counter #to find most common iterations in a list

def buildtracks(filename, dproc,gAli):
 
 dproc.ReadTracksTree(gAli,filename,"nseg>1")
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks

 np = gAli.Npatterns();
 for i in range (np): 
    p = gAli.GetPattern(i);
    ns = p.N();
    for j in range(ns):
      p.GetSegment(j).SetDZ(300);

 return tracks

def trackquality(trackfilepath):
  dproc = ROOT.EdbDataProc()
  gAli = dproc.PVR()
  tracks = buildtracks(trackfilepath, dproc, gAli)
  #list of tracks
  for track in tracks:
    #loop on segments
    nseg = track.N()
    for i in range(nseg):
        segment = track.GetSegment(i)
        mcevent = segment.MCEvt()
        mctrack = segment.MCTrack()

def vertexquality(vertexfilepath):
  inputfile = r.TFile.Open(vertexfilepath,'read')
  inputtree = inputfile.Get("vtx")
  nentries = inputtree.GetEntries()
  for event in inputtree: #loop on the events
    eventlist = event.MCEventID
    trackIDlist = event.MCTrackID
    mothertrackIDlist = event.MCMotherID

    eventcounter = Counter(eventlist)
    mostcommonevent = eventcounter.most_common(1)[0][0]
    eventfrequency = eventcounter.most_common(1)[0][1]
    eventfrequency =  float(eventfrequency)/len(eventlist)

    mothercounter = Counter(mothertrackIDlist)
    mostcommonmother = mothercounter.most_common(1)[0][0]
    motherfrequency = mothercounter.most_common(1)[0][1]
    motherfrequency =  float(motherfrequency)/len(mothertrackIDlist)
    
    hevent.Fill(eventfrequency)
    hmothertrack.Fill(motherfrequency)

hevent = r.TH1F("hevent","Most common Event frequency in vertices",110,0,1.1)
hmothertrack = r.TH1F("hmothertrack","Most common mother track ID frequency in vertices",110,0,1.1)

vertexquality("/home/antonio/Lavoro/Analisi/test_fedra_sim/vertices_MC_modified.root")
cvertexquality = r.TCanvas()
cvertexquality.Divide(1,2)
cvertexquality.cd(1)
hevent.Draw()
cvertexquality.cd(2)
hmothertrack.Draw()