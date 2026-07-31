//Putting the dreaded comparison into plots...(26 Marzo 2026)

void compare_cascade(){
    
    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19;//replace to have multiple years of data taking
    
    TFile *simfile_2018 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2018/nuhistos_Thomas_run1GeVCharm/pythia8_Geant4_charm_1.0_nu.root");    
    TFile *simfile_2026 = TFile::Open("/eos/experiment/ship/user/aiuliano/nuhistos_bkgproductions/bkg2026/makeCascadeTungsten_2026_06_12/all100Runs_Decay_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");
    
    hnu_e_2018->Scale(normship/normsim);
    hnu_e_bar_2018->Scale(normship/normsim);

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2026 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2026 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2026->Scale(normship/normsim);
    hnu_e_bar_2026->Scale(normship/normsim);

    hnu_e_2026->Add(hnu_e_bar_2026);
    hnu_e_2026->SetTitle("electron neutrino from 2026 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Scale(normship/normsim);
    hnu_mu_bar_2018->Scale(normship/normsim);

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2026 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2026 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2026->Scale(normship/normsim);
    hnu_mu_bar_2026->Scale(normship/normsim);

    hnu_mu_2026->Add(hnu_mu_bar_2026);
    hnu_mu_2026->SetTitle("muon neutrino from 2026 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Scale(normship/normsim);
    hnu_tau_bar_2018->Scale(normship/normsim);

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2026 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2026 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2026->Scale(normship/normsim);
    hnu_tau_bar_2026->Scale(normship/normsim);

    hnu_tau_2026->Add(hnu_tau_bar_2026);
    hnu_tau_2026->SetTitle("tau neutrino from 2026 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production from charm cascade");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2026->SetLineColor(kRed);
    hnu_e_2026->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production");
    hnu_e_2026->SetTitle("2026 production");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production from charm cascade");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2026->SetLineColor(kRed);
    hnu_mu_2026->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production");
    hnu_mu_2026->SetTitle("2026 production");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production from charm cascade");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2026->SetLineColor(kRed);
    hnu_tau_2026->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production");
    hnu_tau_2026->SetTitle("2026 production");
    cnutau->BuildLegend();

}

void compare_mbias_nocharm(){

    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19;//replace to have multiple years of data taking

    TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_nocharm/pythia8_Geant4_1.0_c0-19000_nu.root");
    //TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/pythia8_Geant4_1.0_c_nu.root");
    TFile *simfile_2026 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_nocharm/pythia8_Geant4_1.0_c0-157_nu.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");

    hnu_e_2018->Scale(normship/normsim);
    hnu_e_bar_2018->Scale(normship/normsim);

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2026 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2026 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2026->Scale(normship/normsim);
    hnu_e_bar_2026->Scale(normship/normsim);

    hnu_e_2026->Add(hnu_e_bar_2026);
    hnu_e_2026->SetTitle("electron neutrino from 2026 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Scale(normship/normsim);
    hnu_mu_bar_2018->Scale(normship/normsim);

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2026 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2026 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2026->Scale(normship/normsim);
    hnu_mu_bar_2026->Scale(normship/normsim);

    hnu_mu_2026->Add(hnu_mu_bar_2026);
    hnu_mu_2026->SetTitle("muon neutrino from 2026 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Scale(normship/normsim);
    hnu_tau_bar_2018->Scale(normship/normsim);

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2026 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2026 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2026->Scale(normship/normsim);
    hnu_tau_bar_2026->Scale(normship/normsim);

    hnu_tau_2026->Add(hnu_tau_bar_2026);
    hnu_tau_2026->SetTitle("tau neutrino from 2026 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production  (mbias without charm)");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2026->SetLineColor(kRed);
    hnu_e_2026->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production");
    hnu_e_2026->SetTitle("2026 production");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production  (mbias without charm)");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2026->SetLineColor(kRed);
    hnu_mu_2026->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production");
    hnu_mu_2026->SetTitle("2026 production");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production (mbias without charm)");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2026->SetLineColor(kRed);
    hnu_tau_2026->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production");
    hnu_tau_2026->SetTitle("2026 production");
    cnutau->BuildLegend();
}

void compare_mbiasnocharm_pluscascade(){
    
    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19;//replace to have multiple years of data taking
    
    TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_pluscascade_distributions_comparisons/pythia8_Geant4_1.0_withCharm_nu.root");    
    TFile *simfile_2026 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_pluscascade_distributions_comparisons/pythia8_Geant4_1.0_withCharm_nu_makecascade2026_W.root");

    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");
    
    hnu_e_2018->Scale(normship/normsim);
    hnu_e_bar_2018->Scale(normship/normsim);

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    TH1D * hnu_e_2026 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2026 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2026->Scale(normship/normsim);
    hnu_e_bar_2026->Scale(normship/normsim);

    hnu_e_2026->Add(hnu_e_bar_2026);
    hnu_e_2026->SetTitle("electron neutrino from 2026 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Scale(normship/normsim);
    hnu_mu_bar_2018->Scale(normship/normsim);

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2026 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2026 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2026->Scale(normship/normsim);
    hnu_mu_bar_2026->Scale(normship/normsim);

    hnu_mu_2026->Add(hnu_mu_bar_2026);
    hnu_mu_2026->SetTitle("muon neutrino from 2026 production");

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Scale(normship/normsim);
    hnu_tau_bar_2018->Scale(normship/normsim);

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2026 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2026 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2026->Scale(normship/normsim);
    hnu_tau_bar_2026->Scale(normship/normsim);

    hnu_tau_2026->Add(hnu_tau_bar_2026);
    hnu_tau_2026->SetTitle("tau neutrino from 2026 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production from charm cascade");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2026->SetLineColor(kRed);
    hnu_e_2026->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_e_2026->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production from charm cascade");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2026->SetLineColor(kRed);
    hnu_mu_2026->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_mu_2026->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production from charm cascade");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2026->SetLineColor(kRed);
    hnu_tau_2026->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production (mbiasnocharm_pluscascade)");
    hnu_tau_2026->SetTitle("2026 production (mbiasnocharm_pluscascade)");
    cnutau->BuildLegend();

}

void compare_mbiaswithcharm(){

    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19;//replace to have multiple years of data taking
    //double normHanae = 75000* 1923; //NPOT * NFiles
    double normHanae = normsim;

    TFile *simfile_2018 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2018/nuhistos_Thomas_noCharmFalse/pythia8_Geant4_1.0_c0-19000_nu.root");
    //TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/pythia8_Geant4_1.0_c_nu.root");
    //TFile *simfile_2026 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbiaswithcharm/pythia8_Geant4_1.0_c0-157_nu.root");
    //TFile *simfile_2026 = TFile::Open("/eos/experiment/ship/user/htilquin/shared_files/tungsten_check_ship_week_5.root");
    TFile *simfile_2026 = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/nuhistos_bkgproductions/bkg2026/nuhistos_Hanae2026_70GeV_noCharmFalse/pythia8_Geant4_70.0_c0-157_nu.root");
    //nue comparison
    TH1D * hnu_e_2018 = (TH1D*) simfile_2018->Get("1012");
    TH1D * hnu_e_bar_2018 = (TH1D*) simfile_2018->Get("2012");

    hnu_e_2018->Scale(normship/normsim);
    hnu_e_bar_2018->Scale(normship/normsim);

    hnu_e_2018->Add(hnu_e_bar_2018);
    hnu_e_2018->SetTitle("electron neutrino from 2018 production");

    //TH1D * hnu_e_2026 = (TH1D*) simfile_2026->Get("PlaneHA10012");
    TH1D * hnu_e_2026 = (TH1D*) simfile_2026->Get("1012");
    TH1D * hnu_e_bar_2026 = (TH1D*) simfile_2026->Get("2012");

    hnu_e_2026->Scale(normship/normHanae);
    hnu_e_bar_2026->Scale(normship/normHanae);

    hnu_e_2026->Add(hnu_e_bar_2026);
    hnu_e_2026->SetTitle("electron neutrino from 2026 production");

    //numu comparison
    TH1D * hnu_mu_2018 = (TH1D*) simfile_2018->Get("1014");
    TH1D * hnu_mu_bar_2018 = (TH1D*) simfile_2018->Get("2014");

    hnu_mu_2018->Scale(normship/normsim);
    hnu_mu_bar_2018->Scale(normship/normsim);

    hnu_mu_2018->Add(hnu_mu_bar_2018);
    hnu_mu_2018->SetTitle("muon neutrino from 2018 production");

    TH1D * hnu_mu_2026 = (TH1D*) simfile_2026->Get("1014");
    TH1D * hnu_mu_bar_2026 = (TH1D*) simfile_2026->Get("2014");

    hnu_mu_2026->Add(hnu_mu_bar_2026);
    hnu_mu_2026->SetTitle("muon neutrino from 2026 production");

    hnu_mu_2026->Scale(normship/normHanae);
    hnu_mu_bar_2026->Scale(normship/normHanae);

    //nutau comparison
    TH1D * hnu_tau_2018 = (TH1D*) simfile_2018->Get("1016");
    TH1D * hnu_tau_bar_2018 = (TH1D*) simfile_2018->Get("2016");

    hnu_tau_2018->Scale(normship/normsim);
    hnu_tau_bar_2018->Scale(normship/normsim);

    hnu_tau_2018->Add(hnu_tau_bar_2018);
    hnu_tau_2018->SetTitle("tau neutrino from 2018 production");

    TH1D * hnu_tau_2026 = (TH1D*) simfile_2026->Get("1016");
    TH1D * hnu_tau_bar_2026 = (TH1D*) simfile_2026->Get("2016");

    hnu_tau_2026->Scale(normship/normHanae);
    hnu_tau_bar_2026->Scale(normship/normHanae);

    hnu_tau_2026->Add(hnu_tau_bar_2026);
    hnu_tau_2026->SetTitle("tau neutrino from 2026 production");

    gStyle->SetOptStat("nmri");
    //drawing plots
    TCanvas *cnue = new TCanvas("cnue","Electron neutrino momentum at production");
    hnu_e_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_e_2018->Draw("histo");
    hnu_e_2026->SetLineColor(kRed);
    hnu_e_2026->Draw("histo && SAMES");
    cnue->SetLogy();
    
    hnu_e_2018->SetTitle("2018 production (mbias file)");
    hnu_e_2026->SetTitle("2026 production (70 GeV Cut boost)");
    cnue->BuildLegend();

    TCanvas *cnumu = new TCanvas("cnumu","Muon neutrino momentum at production");
    hnu_mu_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_mu_2018->Draw("histo");
    hnu_mu_2026->SetLineColor(kRed);
    hnu_mu_2026->Draw("histo && SAMES");
    cnumu->SetLogy();
    
    hnu_mu_2018->SetTitle("2018 production (mbias file)");
    hnu_mu_2026->SetTitle("2026 production (70 GeV Cut boost)");
    cnumu->BuildLegend();

    TCanvas *cnutau = new TCanvas("cnutau","Tau neutrino momentum at production");
    hnu_tau_2018->GetXaxis()->SetTitle("P [GeV/c]");
    hnu_tau_2018->Draw("histo");
    hnu_tau_2026->SetLineColor(kRed);
    hnu_tau_2026->Draw("histo && SAMES");
    cnutau->SetLogy();
    
    hnu_tau_2018->SetTitle("2018 production (mbias file)");
    hnu_tau_2026->SetTitle("2026 production (70 GeV Cut boost)");
    cnutau->BuildLegend();
}

void compare_mbiasnocharm_pluscascade_2D(){

    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19;//replace to have multiple years of data taking
    
    TFile *simfile_2018 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_pluscascade_distributions_comparisons/pythia8_Geant4_1.0_withCharm_nu.root");    
    TFile *simfile_2026 = TFile::Open("/home/utente/Simulations/nuhistos_comparisons/mbias_pluscascade_distributions_comparisons/pythia8_Geant4_1.0_withCharm_nu_makecascade2026.root");

    TH2D *hnutau2D_2018 = (TH2D*) simfile_2018->Get("1116");
    TH2D *hnutau2D_2026 = (TH2D*) simfile_2026->Get("1116");

    hnutau2D_2018->Scale(normship/normsim);
    hnutau2D_2026->Scale(normship/normsim);

    TCanvas *cnutau2D_2018 = new TCanvas("cnutau2D_2018","ppt_nutau_2018");
    hnutau2D_2018->GetXaxis()->SetTitle("log10p");
    hnutau2D_2018->GetYaxis()->SetTitle("log10pt");
    hnutau2D_2018->Draw("COLZ");

    TCanvas *cnutau2D_2026 = new TCanvas("cnutau2D_2026","ppt_nutau_2026");
    hnutau2D_2026->GetXaxis()->SetTitle("log10p");
    hnutau2D_2026->GetYaxis()->SetTitle("log10pt");
    hnutau2D_2026->Draw("COLZ");



}