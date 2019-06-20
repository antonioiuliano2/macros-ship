import ROOT
import fedrarootlogon


#read tracks with a generic cut
def buildtracks(filename, dproc,gAli, cut = "nseg>1"):
 dproc.ReadTracksTree(gAli,filename,cut)
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks

 np = gAli.Npatterns();
 for i in range (np): 
    p = gAli.GetPattern(i);
    ns = p.N();
    for j in range(ns):
      p.GetSegment(j).SetDZ(300);

 return tracks

#only tracks from a MC event
def buildtrackseventnumber(filename,eventnumber,dproc,gAli, cut = "nseg>1"):

 dproc.ReadTracksTree(gAli,filename,"nseg>1&&s.eMCEvt=={}".format(eventnumber))
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks

 np = gAli.Npatterns();
 for i in range (np): 
    p = gAli.GetPattern(i);
    ns = p.N();
    for j in range(ns):
      p.GetSegment(j).SetDZ(300);

 return tracks
