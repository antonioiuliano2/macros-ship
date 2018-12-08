import ROOT as r
from matplotlib import pyplot as plt

#Study of electron origin in emulsion target for charm cross section measurement (Created at 6 February 2018). Nota: vitale controllare il minimo di BoxPoint.fZ
#f = r.TFile.Open('/home/utente/Simulations/sim_charmdet/pot/Lead6ECCs_gaussian_sigma_0_5/pythia8_Geant4_1_0.001.root')
#f = r.TFile.Open('/home/utente/Simulations/sim_charmdet/charm_events/quickcheckacceptance/ship.conical.Pythia8CharmOnly-TGeant4.root');
f = r.TFile.Open('/home/utente/Simulations/sim_charmdet/pot/Lead6ECCs_nofieldnogaps_G4only_uniform/pythia8_Geant4_1_0.001.root')
tree = f.Get('cbmsim')


nevents = tree.GetEntries();
ninteractions = tree.GetEntries("MCTrack[1].fStartZ < 0"); #number of protons interacted in the target

print nevents, ninteractions

csigmax_all = r.TCanvas()
csigmay_all = r.TCanvas()
csigmax_beam = r.TCanvas()
csigmay_beam = r.TCanvas()
c0 = r.TCanvas()
c1 = r.TCanvas()
def gausbeam():
 csigmax_all.Divide(2,3)
 csigmay_all.Divide(2,3)

 csigmax_all.cd(1)
 #beam at start
 tree.Draw("MCTrack[0].fStartX>>h0x(120,-6,6)")
 h0x = r.gDirectory.Get("h0x")
 h0x.Fit('gaus')
 h0x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(1)
 tree.Draw("MCTrack[0].fStartY>>h0y(100,-5,5)")
 h0y = r.gDirectory.Get("h0y")
 h0y.Fit('gaus')
 h0y.GetXaxis().SetTitle("cm")

 #tracks in emulsion after different ECC
 r.gStyle.SetOptFit(111)
 csigmax_all.cd(2)
 tree.Draw("BoxPoint.fX>>h1x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 20")
 h1x = r.gDirectory.Get("h1x")
 h1x.Fit("gaus")
 h1x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(2)
 tree.Draw("BoxPoint.fY>>h1y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 20")
 h1y = r.gDirectory.Get("h1y")
 h1y.Fit("gaus")
 h1y.GetXaxis().SetTitle("cm")
 
 csigmax_all.cd(3)
 tree.Draw("BoxPoint.fX>>h2x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 40")
 h2x = r.gDirectory.Get("h2x")
 h2x.Fit("gaus")
 h2x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(3)
 tree.Draw("BoxPoint.fY>>h2y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 40")
 h2y = r.gDirectory.Get("h2y")
 h2y.Fit("gaus") 
 h2y.GetYaxis().SetTitle("cm")

 csigmax_all.cd(4)
 tree.Draw("BoxPoint.fX>>h3x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 60")
 h3x = r.gDirectory.Get("h3x")
 h3x.Fit("gaus")
 h3x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(4)
 tree.Draw("BoxPoint.fY>>h3y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 60")
 h3y = r.gDirectory.Get("h3y")
 h3y.Fit("gaus")
 h3y.GetXaxis().SetTitle("cm")

 csigmax_all.cd(5)
 tree.Draw("BoxPoint.fX>>h4x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 80")
 h4x = r.gDirectory.Get("h4x")
 h4x.Fit("gaus")
 h4x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(5)
 tree.Draw("BoxPoint.fY>>h4y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 80")
 h4y = r.gDirectory.Get("h4y")
 h4y.Fit("gaus")
 h4y.GetXaxis().SetTitle("cm")

 csigmax_all.cd(6)
 tree.Draw("BoxPoint.fX>>h5x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 100")
 h5x = r.gDirectory.Get("h5x")
 h5x.Fit("gaus")
 h5x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(6)
 tree.Draw("BoxPoint.fY>>h5y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 100")
 h5y = r.gDirectory.Get("h5y")
 h5y.Fit("gaus")
 h5y.GetXaxis().SetTitle("cm")

def gausonlybeam():
 csigmax_all.Divide(2,3)
 csigmay_all.Divide(2,3)

 csigmax_all.cd(1)
 #beam at start
 tree.Draw("MCTrack[0].fStartX>>h0x(120,-6,6)")
 h0x = r.gDirectory.Get("h0x")
 h0x.Fit('gaus')
 h0x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(1)
 tree.Draw("MCTrack[0].fStartY>>h0y(100,-5,5)")
 h0y = r.gDirectory.Get("h0y")
 h0y.Fit('gaus')
 h0y.GetXaxis().SetTitle("cm")

 #tracks in emulsion after different ECC
 r.gStyle.SetOptFit(111)
 csigmax_all.cd(2)
 tree.Draw("BoxPoint.fX>>h1x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 20 && BoxPoint.fTrackID == 0")
 h1x = r.gDirectory.Get("h1x")
 h1x.Fit("gaus")
 h1x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(2)
 tree.Draw("BoxPoint.fY>>h1y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 20 && BoxPoint.fTrackID == 0")
 h1y = r.gDirectory.Get("h1y")
 h1y.Fit("gaus")
 h1y.GetXaxis().SetTitle("cm")
 
 csigmax_all.cd(3)
 tree.Draw("BoxPoint.fX>>h2x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 40 && BoxPoint.fTrackID == 0")
 h2x = r.gDirectory.Get("h2x")
 h2x.Fit("gaus")
 h2x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(3)
 tree.Draw("BoxPoint.fY>>h2y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 40 && BoxPoint.fTrackID == 0")
 h2y = r.gDirectory.Get("h2y")
 h2y.Fit("gaus") 
 h2y.GetYaxis().SetTitle("cm")

 csigmax_all.cd(4)
 tree.Draw("BoxPoint.fX>>h3x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 60 && BoxPoint.fTrackID == 0")
 h3x = r.gDirectory.Get("h3x")
 h3x.Fit("gaus")
 h3x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(4)
 tree.Draw("BoxPoint.fY>>h3y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 60 && BoxPoint.fTrackID == 0")
 h3y = r.gDirectory.Get("h3y")
 h3y.Fit("gaus")
 h3y.GetXaxis().SetTitle("cm")

 csigmax_all.cd(5)
 tree.Draw("BoxPoint.fX>>h4x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 80 && BoxPoint.fTrackID == 0")
 h4x = r.gDirectory.Get("h4x")
 h4x.Fit("gaus")
 h4x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(5)
 tree.Draw("BoxPoint.fY>>h4y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 80 && BoxPoint.fTrackID == 0")
 h4y = r.gDirectory.Get("h4y")
 h4y.Fit("gaus")
 h4y.GetXaxis().SetTitle("cm")

 csigmax_all.cd(6)
 tree.Draw("BoxPoint.fX>>h5x(120,-6,6)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 100 && BoxPoint.fTrackID == 0")
 h5x = r.gDirectory.Get("h5x")
 h5x.Fit("gaus")
 h5x.GetXaxis().SetTitle("cm")
 csigmay_all.cd(6)
 tree.Draw("BoxPoint.fY>>h5y(100,-5,5)","BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112 && BoxPoint.fDetectorID == 100 && BoxPoint.fTrackID == 0")
 h5y = r.gDirectory.Get("h5y")
 h5y.Fit("gaus")
 h5y.GetXaxis().SetTitle("cm")


def sigmagraph():
 #Writing spreads to plot 
 mm5x = [0.47, 0.660, 0.668,0.730, 0.771, 0.791]
 cm1x = [0.96, 1.01, 1.12, 1.251, 1.308, 1.390]
 cm0_2x = [0.198, 0.324, 0.414, 0.469, 0.537,0.572]
 cm2x = [1.92, 1.97, 2.05, 2.074, 2.116, 2.229]

 mm5y = [0.46, 0.570, 0.682, 0.749, 0.741, 0.813]
 cm0_2y = [0.207, 0.328, 0.404, 0.476, 0.519,0.586] 
 cm1y = [0.97, 1.12, 1.19, 1.235, 1.330, 1.327]
 cm2y = [1.97, 2.04, 2.154, 2.144, 2.216, 2.222]

 #Dividing for sigma_sim
 mm5x = [x /0.5 for x in mm5x] 
 mm5y = [x /0.5 for x in mm5y] 

 cm1x = [x /1.0 for x in cm1x]
 cm1y = [x /1.0 for x in cm1y]

 cm0_2x = [x / 0.2 for x in cm0_2x]
 cm0_2y = [x / 0.2 for x in cm0_2y]
 
 cm2x = [x /2.0 for x in cm2x]
 cm2y = [x /2.0 for x in cm2y] 

 #gstart = r.TGraph(6,index,mm5x)
 plt.plot(mm5x,'m.',label = r'$\sigma_{sim}$ 0.5 cm')
 plt.plot(cm1x,'ko',label = r'$\sigma_{sim}$ 1.0 cm')
 plt.plot(cm0_2x,'y*',label = r'$\sigma_{sim}$ 0.2 cm')
 plt.plot(cm2x,'b.',label = r'$\sigma_{sim}$ 2.0 cm')

 plt.figure(1);
 plt.title('spread of beam at start of each ECC with different sigma')
 plt.xlabel('start of # ECC')
 plt.ylabel('Spread ratio in x')
 plt.legend()

 plt.figure(2);
 plt.plot(mm5y,'m.',label = r'$\sigma_{sim}$ 0.5 cm')
 plt.plot(cm1y,'ko',label = r'$\sigma_{sim}$ 1.0 cm')
 plt.plot(cm0_2y,'y*',label = r'$\sigma_{sim}$ 0.2 cm')
 plt.plot(cm2y,'b*',label = r'$\sigma_{sim}$ 2.0 cm')

 plt.title('spread of beam at start of each ECC with different sigma')
 plt.xlabel('start of # ECC')
 plt.ylabel('Spread ratio in y')
 plt.legend()
 plt.show()

def occupancy():
#use these variables with 'and' does not work as right as it could seem
# noneutralbox = "BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112"
# noneutralspectro = "SpectrometerPoint.fPdgCode != 22 && abs(SpectrometerPoint.fPdgCode) != 2112" 

 #c1 = r.TCanvas()
 c1.Divide(2,2)
 c1.cd(1)
# tree.Draw("SpectrometerPoint.fY:SpectrometerPoint.fX>>h1(4,-2,2,4,-2,2)",noneutralspectro and "SpectrometerPoint.fDetectorID == 101","COLZ")
 tree.Draw("SpectrometerPoint.fY:SpectrometerPoint.fX>>h1(40,-2,2,40,-2,2)","SpectrometerPoint.fPdgCode != 22 && abs(SpectrometerPoint.fPdgCode) != 2112 && SpectrometerPoint.fDetectorID == 111","COLZ")
 c1.cd(2)
 tree.Draw("SpectrometerPoint.fY:SpectrometerPoint.fX>>h2(40,-2,2,40,-2,2)","SpectrometerPoint.fPdgCode != 22 && abs(SpectrometerPoint.fPdgCode) != 2112 && SpectrometerPoint.fDetectorID == 112","COLZ")
 c1.cd(3)
 tree.Draw("SpectrometerPoint.fY:SpectrometerPoint.fX>>h3(40,-2,2,40,-2,2)","SpectrometerPoint.fPdgCode != 22 && abs(SpectrometerPoint.fPdgCode) != 2112 && SpectrometerPoint.fDetectorID == 121","COLZ")
 c1.cd(4)
 tree.Draw("SpectrometerPoint.fY:SpectrometerPoint.fX>>h4(40,-2,2,40,-2,2)","SpectrometerPoint.fPdgCode != 22 && abs(SpectrometerPoint.fPdgCode) != 2112 && SpectrometerPoint.fDetectorID == 122","COLZ")

 c2 = r.TCanvas()
 #MINIMO DI BOXZ PER QUESTA SIMULAZIONE: 39.83
 #tree.Draw("BoxPoint.fZ+39.83>>hbox(10000)", "BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112","histo")
 tree.Draw("BoxPoint.fDetectorID>>hbox", "BoxPoint.fPdgCode != 22 && abs(BoxPoint.fPdgCode) != 2112","histo")

 hbox = r.gDirectory.Get("hbox")
 #hbox.Scale(1./(nevents)) #i clearly separate top and bot in this way
 hbox.Scale(1./(2*nevents))

 hbox.GetXaxis().SetRange(0,hbox.FindBin(20*0.5))
 hbox.Draw()

 print "Maximum in ECC1/2", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))

 hbox.GetXaxis().SetRange(0,hbox.FindBin(20))
 print "Maximum in ECC1", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))
 
 hbox.GetXaxis().SetRange(hbox.FindBin(20),hbox.FindBin(20*2))
 print "Maximum in ECC2", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))
 
 hbox.GetXaxis().SetRange(hbox.FindBin(20*2),hbox.FindBin(20*3))
 print "Maximum in ECC3", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))

 hbox.GetXaxis().SetRange(hbox.FindBin(20*3),hbox.FindBin(20*4))
 print "Maximum in ECC4", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))

 hbox.GetXaxis().SetRange(hbox.FindBin(20*4),hbox.FindBin(20*5))
 print "Maximum in ECC5", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))

 hbox.GetXaxis().SetRange(hbox.FindBin(20*5),hbox.FindBin(20*6))
 print "Maximum in ECC6", hbox.GetMaximum(), "Numero di protoni da esporre: ", 1000./(hbox.GetMaximum()/(125*99))

 hbox.GetXaxis().SetRange(0,hbox.FindBin(20*6))
# hbox.GetXaxis().SetRange(0,122)
 
 c0.cd(1)
 hbox.Draw()
 l0 = r.TLine(10*1,0,10*1,1000);
 l1 = r.TLine(20*1,0,20*1,1000);
 l2 = r.TLine(20*2,0,20*2,1000);
 l3 = r.TLine(20*3,0,20*3,1000);
 l4 = r.TLine(20*4,0,20*4,1000);
 l5 = r.TLine(20*5,0,20*5,1000);
 l6 = r.TLine(20*6,0,20*6,1000);

 l0.Draw("SAME")
 l1.Draw("SAME")
 l2.Draw("SAME")
 l3.Draw("SAME")
 l4.Draw("SAME")
 l5.Draw("SAME")
 l6.Draw("SAME")

 hboxcut = hbox.Clone('hboxcut') #highlighting selection
 hboxcut.SetFillColor(r.kYellow)
 image = r.TImage.Create() #image to be extracted from canvas and saved each time

 hboxcut.GetXaxis().SetRange(0,hboxcut.FindBin(10)-1)
 hboxcut.Draw("SAME&&HISTO")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC0.png");
 hboxcut.GetXaxis().SetRange(10,hboxcut.FindBin(20)-1)
 #hboxcut.Draw("SAME")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC1.png");
 hboxcut.GetXaxis().SetRange(hboxcut.FindBin(20*1),hboxcut.FindBin(20*2)-1)
 #hboxcut.Draw("SAME")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC2.png");
 hboxcut.GetXaxis().SetRange(hboxcut.FindBin(20*2),hboxcut.FindBin(20*3)-1)
 #hboxcut.Draw("SAME")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC3.png");
 hboxcut.GetXaxis().SetRange(hboxcut.FindBin(20*3),hboxcut.FindBin(20*4)-1)
 #hboxcut.Draw("SAME")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC4.png");
 hboxcut.GetXaxis().SetRange(hboxcut.FindBin(20*4)+1,hboxcut.FindBin(20*5))
 #hboxcut.Draw("SAME")
 c2.Draw()
 c2.Update()
 image.FromPad(c2);  
 image.WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/OccupancyECC5.png");

 h1 = r.gDirectory.Get("h1")
 h2 = r.gDirectory.Get("h2")
 h3 = r.gDirectory.Get("h3")
 h4 = r.gDirectory.Get("h4")

 #h1.Scale(1./ninteractions)
 #h2.Scale(1./ninteractions)
 #h3.Scale(1./ninteractions)
 #h4.Scale(1./ninteractions)
 print h2.GetEntries()
 #r.gApplication.Run() #leaving canvas open so I can work with gui





def studyorigin():
 #defining histograms
 hmumpdg = r.TH1I("hmumpdg","Pdg Codes mother", 100, -50, 50)
 htrackmomentum = r.TH1D("htrackmomentum","Momentum of electrons in detectorid 22", 100, 0, 0.1)
 hmummomentum = r.TH1D("hmummomentum","Momentum mother", 100, 0, 0.1)
 hstartz = r.TH1D("hstartZ", "Track start z", 1000, -90,-40)
 hstartzvsmomentum = r.TH2D("hstartzvsmomentum", "Startz vs momentum", 1000, -90, -40, 100, 0, 0.1)
 
 for i in range(100): #loop on events
  if i%10 == 0: print i
  tree.GetEntry(i)
  tracks = tree.MCTrack
  boxpoints = tree.BoxPoint
 # for track in tracks:

  for boxpoint in boxpoints: #loop on tracks
   if boxpoint.GetDetectorID() == 22 and boxpoint.GetTrackID() > 0 and abs(boxpoint.PdgCode() == 11):
    #print boxpoint.GetTrackID(), boxpoint.PdgCode():
    trackid = boxpoint.GetTrackID()
    track = tracks[trackid]
    mumid = track.GetMotherId()
    if mumid > 0:
     htrackmomentum.Fill(track.GetP())
     hstartz.Fill(track.GetStartZ())
     hstartzvsmomentum.Fill(track.GetStartZ(), track.GetP())
     mumtrack = tracks[mumid]
     hmumpdg.Fill(mumtrack.GetPdgCode())
     hmummomentum.Fill(mumtrack.GetP())
   

 c1 = r.TCanvas()
 hmumpdg.Draw() 

 c2 = r.TCanvas()
 hmumpdg.GetXaxis().SetTitle("GeV/c")
 hmummomentum.Draw()

 c3 = r.TCanvas()
 hstartz.GetXaxis().SetTitle("cm")
 hstartz.Draw()

 c4 = r.TCanvas()
 htrackmomentum.GetXaxis().SetTitle("GeV/c")
 htrackmomentum.Draw()

 c5 = r.TCanvas()
 hstartzvsmomentum.GetXaxis().SetTitle("cm")
 hstartzvsmomentum.GetYaxis().SetTitle("GeV/c")
 hstartzvsmomentum.Draw()



