/*# Comparing neutrino spectra with and without cascade

First, we access the 2018 backgorund files from EOS. From the README:

histograms with neutrino momentum:

 * pythia8_Geant4_10.0_c_nu.root (mbias without charm
 
 * pythia8_Geant4_charm_nu_10.0.root (charm from cascade
 
 * pythia8_Geant4_10.0_withCharm_nu.root  (sum of above

 * pythia8_Geant4_1.0_c_nu.root (mbias without charm
 
 * pythia8_Geant4_charm_nu_1.0.root (charm from cascade
 
 * pythia8_Geant4_1.0_withCharm_nu.root  (sum of above
 
We consider here the cut at **1 GeV** */


void cascadecomparison(){

 // 1D P histograms

 TFile *neutrino_withcharm = TFile::Open("pythia8_Geant4_1.0_withCharm_nu.root");
 TFile *neutrino_onlycharm = TFile::Open("pythia8_Geant4_charm_nu_1.0.root");
 TFile *neutrino_nocharm = TFile::Open("pythia8_Geant4_1.0_c_nu.root");

 TH1D *hpnutau_withcharm = (TH1D*) neutrino_withcharm->Get("1016");
 TH1D *hpnutaubar_withcharm = (TH1D*) neutrino_withcharm->Get("2016");
 
 TH1D *hpnutau_nocharm = (TH1D*) neutrino_nocharm->Get("16");
 TH1D *hpnutaubar_nocharm = (TH1D*) neutrino_nocharm->Get("-16");

 TH1D *hpnutau_onlycharm = (TH1D*) neutrino_onlycharm->Get("16");
 TH1D *hpnutaubar_onlycharm = (TH1D*) neutrino_onlycharm->Get("-16");

 //hpnutau_onlycharm->Add(hpnutaubar_onlycharm);
 //hpnutau_withcharm->Add(hpnutaubar_withcharm);
 //hpnutau_nocharm->Add(hpnutaubar_nocharm);

 TCanvas *cP = new TCanvas();
 hpnutau_withcharm->Draw();
 hpnutau_withcharm->SetTitle("Sum of above;P[GeV/c]");
 hpnutau_nocharm->SetLineColor(kRed);
 hpnutau_nocharm->SetTitle("mbias without charm");
 hpnutau_onlycharm->SetLineColor(kGreen);
 hpnutau_onlycharm->SetTitle("charm from cascade");
 hpnutau_onlycharm->Draw("SAMES");
 hpnutau_nocharm->Draw("SAMES");
 cP->SetLogy();
 cP->BuildLegend();

 hpnutau_withcharm->SetTitle("Nu tau momentum;P[GeV/c]");

 //passing at 2D PPT histograms
}
