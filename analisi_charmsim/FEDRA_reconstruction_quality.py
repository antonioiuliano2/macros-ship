#test tracking and reconstruction with MC information

import ROOT as r
import fedrarootlogon
#from collections import Counter #to find most common iterations in a list

def countX(lst, x): 
    count = 0
    for ele in lst: 
        if (ele == x): 
            count = count + 1
    return count 

def most_frequent(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = countX(List, i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num,counter 

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
  dproc = r.EdbDataProc()
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

    mostcommonevent = most_frequent(eventlist)[0]
    eventfrequency = most_frequent(eventlist)[1]
    eventfrequency =  float(eventfrequency)/len(eventlist)

#    mothercounter = (mothertrackIDlist)
    mostcommonmother = most_frequent(mothertrackIDlist)[0]
    motherfrequency = most_frequent(mothertrackIDlist)[1]
    motherfrequency =  float(motherfrequency)/len(mothertrackIDlist)
    
    hevent.Fill(eventfrequency)
    hmothertrack.Fill(motherfrequency)

hevent = r.TH1F("hevent","Most common Event frequency in vertices",110,0,1.1)
hmothertrack = r.TH1F("hmothertrack","Most common mother track ID frequency in vertices",110,0,1.1)

vertexquality("verticesandtracks_2500events.root")
cvertexquality = r.TCanvas()
cvertexquality.Divide(1,2)
cvertexquality.cd(1)
hevent.Draw()
cvertexquality.cd(2)
hmothertrack.Draw()
