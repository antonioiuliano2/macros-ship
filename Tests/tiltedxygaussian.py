'''generating xy according to 2 gaussians, then applying angular transformation'''
import ROOT as r

sigmax = 0.137
sigmay = 0.5

#theta = + r.TMath.Pi()/6.
theta = + 0.20
npoints = int(1E+5)

hxy = r.TH2D("hxy","Produced positions;x[cm];y[cm]",30,-1.5,1.5,30,-1.5,1.5)
hx1y1 = r.TH2D("hx1y1","Tilted positions;x[cm];y[cm]",30,-1.5,1.5,30,-1.5,1.5)
hdt = r.TH1D("hdt","DeltaTimeStamp",200,0,200)
ht = r.TH1D("ht","TimeStamp",200,0,200000)
randomgen = r.TRandom3()
for ipoint in range(npoints):
    #generating initial pair
    x = randomgen.Gaus(0,sigmax)
    y = randomgen.Gaus(0,sigmay)
    #applying rotation
    vec = r.TVector3(x,y,0.)
    vec.RotateZ(theta)
    #x1 = x* r.TMath.Cos(theta) - y *r.TMath.Sin(theta)
    #y1 =  x * r.TMath.Sin(theta) + y*r.TMath.Cos(theta)
    x1 = vec[0]
    y1 = vec[1]
    hxy.Fill(x,y) 
    hx1y1.Fill(x1,y1)
t = 0.
ngauspoints = 15000
for ipoint in range(ngauspoints):
 dt = r.gRandom.Exp(7.721)
 t = t+dt
 hdt.Fill(dt)
 ht.Fill(t)
 
c = r.TCanvas()
c.Divide(1,2)
c.cd(1)
hxy.Draw("COLZ")
c.cd(2)
hx1y1.Draw("COLZ")

ct = r.TCanvas()
ct.Divide(1,2)
ct.cd(1)
ht.Draw()
ct.cd(2)
hdt.Draw()