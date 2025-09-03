//obtain the statistics error from the number of neutrinos, compute errors on cross section (A Iuliano, 3 September 2025)

void statistics_error_xsec(){
 //from getnuestimation(), the part where I get the histograms from file
 TString prepath("/home/utente/Simulations/nuyield_shipecn3/2025_08_28_nuyield_SND_EmuTargetSiliconTarget/"); //launch from here, at least for now

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 4e+19; //replace to have multiple years of data taking
 //double normship = 6e+20; //15 years
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
 numuspectrum->Scale(normship/normsim);
 nutauspectrum->Scale(normship/normsim);

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
 numubarspectrum->Scale(normship/normsim);
 nutaubarspectrum->Scale(normship/normsim);

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

 //then, we retrieve the cross section
 const Int_t A = 184; //Tungsten as passive material
 const char* nu="nu_mu";
 const char* intmode="dis_cc";
 TFile *xsec = TFile::Open("nu_xsec_TungstenSHIP_emax400.root");
 //N.D.R! These are for 10^-38 cm2, and they are for all protons and neutrons. We need to divide it x 184 to get it per Nucleon
 TGraph *hxsec_p = (TGraph*) xsec->Get(Form("%s_W184/%s_p",nu,intmode));   
 TGraph *hxsec_n = (TGraph*) xsec->Get(Form("%s_W184/%s_n",nu,intmode));

 TGraph *hxsec_bar_p = (TGraph*) xsec->Get(Form("%s_bar_W184/%s_p",nu,intmode));   
 TGraph *hxsec_bar_n = (TGraph*) xsec->Get(Form("%s_bar_W184/%s_n",nu,intmode));
 //we also need to divide for energy
 TGraphErrors *hxsec_N_overE = new TGraphErrors();
 TGraphErrors *hxsec_bar_N_overE = new TGraphErrors();

 double minenergy = 10.; //cross section is not well defined under 4 GeV
 double maxenergy = 200.; //not much point for statistics error
 int ipoint=0;
 for (int ibin=1; ibin<=nbins;ibin++){
  double energy = numuspectrum->GetXaxis()->GetBinCenter(ibin);
  if (energy<minenergy || energy>maxenergy) continue;
  double xsec = (hxsec_p->Eval(energy) + hxsec_n->Eval(energy))/(A * energy);
  double delta_xsec = numuspectrum->GetBinError(ibin)/numuspectrum->GetBinContent(ibin)*xsec; 

  hxsec_N_overE->AddPoint(energy,xsec);
  hxsec_N_overE->SetPointError(ipoint,numuspectrum->GetBinWidth(ibin)/2.,delta_xsec);

  double xsec_bar = (hxsec_bar_p->Eval(energy) + hxsec_bar_n->Eval(energy))/(A * energy);
  double delta_xsec_bar = numubarspectrum->GetBinError(ibin)/numubarspectrum->GetBinContent(ibin)*xsec_bar; 

  hxsec_bar_N_overE->AddPoint(energy,xsec_bar);
  hxsec_bar_N_overE->SetPointError(ipoint,numubarspectrum->GetBinWidth(ibin)/2.,delta_xsec_bar);

  ipoint++;
 }

 TCanvas *cxsec = new TCanvas();
 hxsec_N_overE->Draw();
 hxsec_N_overE->GetYaxis()->SetRangeUser(0.1,0.9);
 hxsec_bar_N_overE->Draw("SAME");

 cxsec->SetLogx();

}