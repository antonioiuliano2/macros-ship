//obtain the statistics error from the number of neutrinos, compute errors on cross section (A Iuliano, 3 September 2025)
TCanvas * computexsec_errors(TFile* xsec_file, const char* nu, const char* intmode, TH1D *nutauspectrum, TH1D *nutaubarspectrum, const char* charmmode=""){
 TGraphErrors *hxsec_N_overE = new TGraphErrors();
 TGraphErrors *hxsec_bar_N_overE = new TGraphErrors();
 //N.D.R! These are for 10^-38 cm2, and they are for all protons and neutrons. We need to divide it x 184 to get it per Nucleon
 const Int_t A = 184; //Tungsten as passive material
 TGraph *hxsec_p = (TGraph*) xsec_file->Get(Form("%s_W184/%s_p%s",nu,intmode,charmmode));   
 TGraph *hxsec_n = (TGraph*) xsec_file->Get(Form("%s_W184/%s_n%s",nu,intmode,charmmode));

 TGraph *hxsec_bar_p = (TGraph*) xsec_file->Get(Form("%s_bar_W184/%s_p%s",nu,intmode,charmmode));   
 TGraph *hxsec_bar_n = (TGraph*) xsec_file->Get(Form("%s_bar_W184/%s_n%s",nu,intmode,charmmode));
 //we also need to divide for energy

 double minenergy = 10.; //cross section is not well defined under 4 GeV
 double maxenergy = 200.; //not much point for statistics error
 int ipoint=0;
 int nbins = nutauspectrum->GetNbinsX();
 for (int ibin=1; ibin<=nbins;ibin++){
  double energy = nutauspectrum->GetXaxis()->GetBinCenter(ibin);
  if (energy<minenergy || energy>maxenergy) continue;
  double xsec = (hxsec_p->Eval(energy) + hxsec_n->Eval(energy))/(A * energy);
  double delta_xsec = nutauspectrum->GetBinError(ibin)/nutauspectrum->GetBinContent(ibin)*xsec; 

  hxsec_N_overE->AddPoint(energy,xsec);
  hxsec_N_overE->SetPointError(ipoint,nutauspectrum->GetBinWidth(ibin)/2.,delta_xsec);

  double xsec_bar = (hxsec_bar_p->Eval(energy) + hxsec_bar_n->Eval(energy))/(A * energy);
  double delta_xsec_bar = nutaubarspectrum->GetBinError(ibin)/nutaubarspectrum->GetBinContent(ibin)*xsec_bar; 

  hxsec_bar_N_overE->AddPoint(energy,xsec_bar);
  hxsec_bar_N_overE->SetPointError(ipoint,nutaubarspectrum->GetBinWidth(ibin)/2.,delta_xsec_bar);

  ipoint++;
 }

 TCanvas *cxsec = new TCanvas(Form("c_%s_%s%s",nu,intmode,charmmode),Form("%s_%s%s",nu,intmode,charmmode));
 hxsec_N_overE->GetYaxis()->SetRangeUser(0.1,0.9);
 hxsec_N_overE->Draw();
 hxsec_bar_N_overE->SetLineColor(kBlue);
 hxsec_bar_N_overE->Draw("SAME");
 hxsec_N_overE->SetTitle("#nu;E[GeV];#sigma_{#nuN}(E)/E[10^{-38}cm^{2}/GeV]");
 hxsec_bar_N_overE->SetTitle("anti-#nu;E[GeV];#sigma_{#nuN}(E)/E[10^{-38}cm^{2}/GeV]");

 cxsec->SetLogx();
 cxsec->BuildLegend();
 return cxsec;
}

void statistics_error_xsec(){
 //from getnuestimation(), the part where I get the histograms from file
 TString prepath("/home/utente/Simulations/nuyield_shipecn3/2025_08_28_nuyield_SND_EmuTargetSiliconTarget/"); //launch from here, at least for now

 double numu_deteff = 0.54; //efficiency old SHIP@ECN3 document study, 2023_07_18
 double nutau_deteff = 0.30; // only in the muon channel, need also to consider the branching ratio
 double tau_inmuon_br = 0.174; // don't forget, I am with you in the tau

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 //double normship = 4e+19; //replace to have multiple years of data taking
 double normship = 6e+20; //15 years
 TFile *nuefile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_e_dis_cc.root").Data());
 TFile *numufile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_mu_dis_cc.root").Data());
 TFile *nutaufile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_tau_dis_cc.root").Data());
 TFile *nuebarfile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_e_bar_dis_cc.root").Data());
 TFile *numubarfile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_mu_bar_dis_cc.root").Data());
 TFile *nutaubarfile = TFile::Open(prepath+TString("logbinning_plots_2/results_nu_tau_bar_dis_cc.root").Data());


 TH1D *nuespectrum = (TH1D*) nuefile->Get("hspectrum_nu_e_intdis_cc");
 TH1D *numuspectrum = (TH1D*) numufile->Get("hspectrum_nu_mu_intdis_cc");
 TH1D *nutauspectrum = (TH1D*) nutaufile->Get("hspectrum_nu_tau_intdis_cc");

 nuespectrum->Scale(normship/normsim);
 numuspectrum->Scale(normship/normsim * numu_deteff);
 nutauspectrum->Scale(normship/normsim * nutau_deteff * tau_inmuon_br);

 //printing yields to screen
 cout<<"nue: Mean energy:\t"<<nuespectrum->GetMean()<<" N:\t"<<nuespectrum->Integral()<<endl;
 cout<<endl;
 cout<<"numu: Mean energy:\t"<<numuspectrum->GetMean()<<" N:\t"<<numuspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"nutau: Mean energy:\t"<<nutauspectrum->GetMean()<<" N:\t"<<nutauspectrum->Integral()<<endl;
 cout<<endl;

 TH1D *nuebarspectrum = (TH1D*) nuebarfile->Get("hspectrum_nu_e_bar_intdis_cc");
 TH1D *numubarspectrum = (TH1D*) numubarfile->Get("hspectrum_nu_mu_bar_intdis_cc");
 TH1D *nutaubarspectrum = (TH1D*) nutaubarfile->Get("hspectrum_nu_tau_bar_intdis_cc");

 nuebarspectrum->Scale(normship/normsim);
 numubarspectrum->Scale(normship/normsim * numu_deteff);
 nutaubarspectrum->Scale(normship/normsim * nutau_deteff * tau_inmuon_br);

 //printing yields to screen
 cout<<"nuebar: Mean energy:\t"<<nuebarspectrum->GetMean()<<" N:\t"<<nuebarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"numubar: Mean energy:\t"<<numubarspectrum->GetMean()<<" N:\t"<<numubarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"nutaubar: Mean energy:\t"<<nutaubarspectrum->GetMean()<<" N:\t"<<nutaubarspectrum->Integral()<<endl;
 cout<<endl;

 int nbins = nutauspectrum->GetNbinsX();
 for (int ibin=1; ibin<=nbins;ibin++){
  nuespectrum->SetBinError(ibin,TMath::Sqrt(nuespectrum->GetBinContent(ibin)));
  nuebarspectrum->SetBinError(ibin,TMath::Sqrt(nuebarspectrum->GetBinContent(ibin)));

  numuspectrum->SetBinError(ibin,TMath::Sqrt(numuspectrum->GetBinContent(ibin)));
  numubarspectrum->SetBinError(ibin,TMath::Sqrt(numubarspectrum->GetBinContent(ibin)));

  nutauspectrum->SetBinError(ibin,TMath::Sqrt(nutauspectrum->GetBinContent(ibin)));
  nutaubarspectrum->SetBinError(ibin,TMath::Sqrt(nutaubarspectrum->GetBinContent(ibin)));
 }

 TCanvas *cnu = new TCanvas();
 nutauspectrum->SetLineColor(kGreen);
 nutauspectrum->Draw();
 nutaubarspectrum->Draw("SAMES");

 nutauspectrum->SetTitle("nu_tau");
 nutaubarspectrum->SetTitle("anti-nu_tau");
 cnu->BuildLegend();

 TFile *xsec_file = TFile::Open((prepath+TString("nu_xsec_TungstenSHIP_emax400.root")).Data());
 //then, we retrieve the cross section
 const char* nu="nu_tau";
 const char* intmode="dis_cc";
  
 TCanvas *cxsec_nutau = computexsec_errors(xsec_file, "nu_tau", "dis_cc", nutauspectrum, nutaubarspectrum);
 TCanvas *cxsec_numu = computexsec_errors(xsec_file, "nu_mu", "dis_cc", numuspectrum, numubarspectrum);
 TCanvas *cxsec_numu_charm = computexsec_errors(xsec_file, "nu_mu", "dis_cc", numuspectrum, numubarspectrum,"_charm");

}


