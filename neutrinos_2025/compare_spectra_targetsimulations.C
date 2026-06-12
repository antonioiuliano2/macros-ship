//Putting the dreaded comparison into plots...(26 Marzo 2026)

void compare_cascade(){
    
    
    TFile *simfile_2018 = TFile::Open("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_charm_nu_1.0.root");    
    TFile *simfile_2026 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2026/nuhistos_Kirill2026_CharmCascade/pythia8_Geant4_charm_1.0_nu.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2025 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2025 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2025->Add(hnu_e_bar_2025);
    hnu_e_2025->SetTitle("electron neutrino from 2025 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2025 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2025 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2025->Add(hnu_mu_bar_2025);
    hnu_mu_2025->SetTitle("muon neutrino from 2025 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2025 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2025 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2025->Add(hnu_tau_bar_2025);
    hnu_tau_2025->SetTitle("tau neutrino from 2025 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production from charm cascade");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2025->SetLineColor(kRed);
    hnu_e_2025->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production");
    hnu_e_2025->SetTitle("2026 production");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production from charm cascade");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2025->SetLineColor(kRed);
    hnu_mu_2025->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production");
    hnu_mu_2025->SetTitle("2026 production");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production from charm cascade");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2025->SetLineColor(kRed);
    hnu_tau_2025->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production");
    hnu_tau_2025->SetTitle("2026 production");
    cnutau->BuildLegend();

}

void compare_mbias_nocharm(){

    TFile *simfile_2018 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2018/nuhistos_Thomas_noCharmTrue/pythia8_Geant4_1.0_c0-19000_nu.root");
    //TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/pythia8_Geant4_1.0_c_nu.root");
    TFile *simfile_2026 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2026/nuhistos_Hanae2026_noCharmTrue/pythia8_Geant4_1.0_c0-157_nu.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2025 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2025 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2025->Add(hnu_e_bar_2025);
    hnu_e_2025->SetTitle("electron neutrino from 2025 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2025 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2025 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2025->Add(hnu_mu_bar_2025);
    hnu_mu_2025->SetTitle("muon neutrino from 2025 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2025 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2025 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2025->Add(hnu_tau_bar_2025);
    hnu_tau_2025->SetTitle("tau neutrino from 2025 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production  (mbias without charm)");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2025->SetLineColor(kRed);
    hnu_e_2025->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production");
    hnu_e_2025->SetTitle("2026 production");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production  (mbias without charm)");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2025->SetLineColor(kRed);
    hnu_mu_2025->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production");
    hnu_mu_2025->SetTitle("2026 production");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production (mbias without charm)");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2025->SetLineColor(kRed);
    hnu_tau_2025->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production");
    hnu_tau_2025->SetTitle("2026 production");
    cnutau->BuildLegend();
}

void compare_mbiasnocharm_pluscascade(){
    
    
    TFile *simfile_2018 = TFile::Open("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_1.0_withCharm_nu.root");    
    TFile *simfile_2026 = TFile::Open("/eos/experiment/ship/user/aiuliano/nuhistos_bkgproductions/bkg2026/nuhistos_Hanae_Kirill_2026_merged/pythia8_Geant4_1.0_withCharm_nu.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2025 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2025 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2025->Add(hnu_e_bar_2025);
    hnu_e_2025->SetTitle("electron neutrino from 2025 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2025 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2025 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2025->Add(hnu_mu_bar_2025);
    hnu_mu_2025->SetTitle("muon neutrino from 2025 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2025 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2025 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2025->Add(hnu_tau_bar_2025);
    hnu_tau_2025->SetTitle("tau neutrino from 2025 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production from charm cascade");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2025->SetLineColor(kRed);
    hnu_e_2025->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_e_2025->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production from charm cascade");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2025->SetLineColor(kRed);
    hnu_mu_2025->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_mu_2025->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production from charm cascade");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2025->SetLineColor(kRed);
    hnu_tau_2025->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_tau_2025->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnutau->BuildLegend();

}

void compare_spectra_targetsimulations(){

    TFile *simfile_2018 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2018/nuhistos_Thomas_noCharmFalse_checkmother/pythia8_Geant4_1.0_c0-19000_nu.root");
    //TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/pythia8_Geant4_1.0_c_nu.root");
    TFile *simfile_2026 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2026/nuhistos_Hanae2026_noCharmFalse_checkmother/pythia8_Geant4_1.0_c0-157_nu.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("41012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("42012");

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2025 = (TH1D*) simfile_2026->Get("41012");
    TH1D * hnu_e_bar_2025 = (TH1D*) simfile_2026->Get("42012");

    hnu_e_2025->Add(hnu_e_bar_2025);
    hnu_e_2025->SetTitle("electron neutrino from 2025 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("41014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("42014");

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2025 = (TH1D*) simfile_2026->Get("41014");
    TH1D * hnu_mu_bar_2025 = (TH1D*) simfile_2026->Get("42014");

    hnu_mu_2025->Add(hnu_mu_bar_2025);
    hnu_mu_2025->SetTitle("muon neutrino from 2025 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("41016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("42016");

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2025 = (TH1D*) simfile_2026->Get("41016");
    TH1D * hnu_tau_bar_2025 = (TH1D*) simfile_2026->Get("42016");

    hnu_tau_2025->Add(hnu_tau_bar_2025);
    hnu_tau_2025->SetTitle("tau neutrino from 2025 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2025->SetLineColor(kRed);
    hnu_e_2025->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production (mbias file)");
    hnu_e_2025->SetTitle("2026 production");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2025->SetLineColor(kRed);
    hnu_mu_2025->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production (mbias file)");
    hnu_mu_2025->SetTitle("2026 production");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2025->SetLineColor(kRed);
    hnu_tau_2025->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production (mbias file)");
    hnu_tau_2025->SetTitle("2026 production");
    cnutau->BuildLegend();
}