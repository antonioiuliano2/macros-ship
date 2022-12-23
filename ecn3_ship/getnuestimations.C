//just draw histograms together, with colors as in SHiP CDS plots, get them together
void getnuresults(){

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 2e+20; //reference for five years of SND DataTaking

 TFile *nuefile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_e_dis_cc.root");
 TFile *numufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_mu_dis_cc.root");
 TFile *nutaufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_tau_dis_cc.root");
 TFile *nuebarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_e_bar_dis_cc.root");
 TFile *numubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_mu_bar_dis_cc.root");
 TFile *nutaubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/plots/results_nu_tau_bar_dis_cc.root");


 TH1D *nuespectrum = (TH1D*) nuefile->Get("hspectrum_nu_e_intdis_cc");
 TH1D *numuspectrum = (TH1D*) numufile->Get("hspectrum_nu_mu_intdis_cc");
 TH1D *nutauspectrum = (TH1D*) nutaufile->Get("hspectrum_nu_tau_intdis_cc");

 nuespectrum->Scale(normship/normsim);
 numuspectrum->Scale(normship/normsim);
 nutauspectrum->Scale(normship/normsim);

//printing yields to screen
 cout<<"N nue"<<nuespectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N numu"<<numuspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N nutau"<<nutauspectrum->Integral()<<endl;
 cout<<endl;

 TH1D *nuebarspectrum = (TH1D*) nuebarfile->Get("hspectrum_nu_e_bar_intdis_cc");
 TH1D *numubarspectrum = (TH1D*) numubarfile->Get("hspectrum_nu_mu_bar_intdis_cc");
 TH1D *nutaubarspectrum = (TH1D*) nutaubarfile->Get("hspectrum_nu_tau_bar_intdis_cc");

 nuebarspectrum->Scale(normship/normsim);
 numubarspectrum->Scale(normship/normsim);
 nutaubarspectrum->Scale(normship/normsim);

 //printing yields to screen
 cout<<"N nuebar"<<nuebarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N numubar"<<numubarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N nutaubar"<<nutaubarspectrum->Integral()<<endl;
 cout<<endl;

 //combining nu and antinu together
 nuespectrum->Add(nuebarspectrum);
 numuspectrum->Add(numubarspectrum);
 nutauspectrum->Add(nutaubarspectrum);

 TCanvas *cnu = new TCanvas();
 numuspectrum->SetLineColor(kRed);
 numuspectrum->Draw("hist");
 nuespectrum->SetLineColor(kBlue);
 nuespectrum->Draw("hist&&SAMES");
 nutauspectrum->SetLineColor(kGreen);
 nutauspectrum->Draw("hist&&SAMES");


}

void getnuincomingresults(){
 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 2e+20; //reference for five years of SND DataTaking

 TFile *nufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/neutrinos_detector.root");

 TH1D *nuespectrum = (TH1D*) nufile->Get("hnu_e");
 TH1D *numuspectrum = (TH1D*) nufile->Get("hnu_mu");
 TH1D *nutauspectrum = (TH1D*) nufile->Get("hnu_tau");

 nuespectrum->Scale(normship/normsim);
 numuspectrum->Scale(normship/normsim);
 nutauspectrum->Scale(normship/normsim);
 //printing yields to screen
 cout<<"N nue"<<nuespectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N numu"<<numuspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N nutau"<<nutauspectrum->Integral()<<endl;
 cout<<endl;

 TH1D *nuebarspectrum = (TH1D*) nufile->Get("hnu_e_bar");
 TH1D *numubarspectrum = (TH1D*) nufile->Get("hnu_mu_bar");
 TH1D *nutaubarspectrum = (TH1D*) nufile->Get("hnu_tau_bar");

 nuebarspectrum->Scale(normship/normsim);
 numubarspectrum->Scale(normship/normsim);
 nutaubarspectrum->Scale(normship/normsim);


 //printing yields to screen
 cout<<"N nuebar"<<nuebarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N numubar"<<numubarspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"N nutaubar"<<nutaubarspectrum->Integral()<<endl;
 cout<<endl;

 //combining nu and antinu together
 nuespectrum->Add(nuebarspectrum);
 numuspectrum->Add(numubarspectrum);
 nutauspectrum->Add(nutaubarspectrum);

 TCanvas *cnu = new TCanvas();
 numuspectrum->SetLineColor(kRed);
 numuspectrum->Draw("hist");
 nuespectrum->SetLineColor(kBlue);
 nuespectrum->Draw("hist&&SAMES");
 nutauspectrum->SetLineColor(kGreen);
 nutauspectrum->Draw("hist&&SAMES");

}