import ROOT as r
import fedrarootlogon

dproc = r.EdbDataProc()
gAli = dproc.PVR()

def buildtracks(filename,cut):
 print cut
 dproc.ReadTracksTree(gAli,filename,cut)
 gAli.FillCell(30,30,0.009,0.009)
 tracks = gAli.eTracks
 
 return tracks

xmin = 20000.
ymin = 20000.
xmax = 90000.
ymax = 90000.

nnumbers = 100

for i in range(nnumbers):
 #generating a position in the space
 vx = r.gRandom.Uniform(xmin,xmax)
 vy = r.gRandom.Uniform(ymin,ymax)

 print "generated point in ", vx, vy
 #reading tracks in a 6 mm x 6 mm region around this point
 buildtracks("/ship/CHARM2018/CH1-R6/b000001/b000001.0.0.0.trk.root",("nseg>1&&TMath::Abs(s[0].eX-%d)<6000&&TMath::Abs(s[0].eY-%d)<6000"%(vx,vy)))

 gAli.ResetTracks()
