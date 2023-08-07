//just draw histograms together, with colors as in SHiP CDS plots, get them together
void getnuresults(){

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship =6e+20; //replace to have multiple years of data taking

 TFile *nuefile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_e_dis_cc.root");
 TFile *numufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_mu_dis_cc.root");
 TFile *nutaufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_tau_dis_cc.root");
 TFile *nuebarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_e_bar_dis_cc.root");
 TFile *numubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_mu_bar_dis_cc.root");
 TFile *nutaubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_tau_bar_dis_cc.root");


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

void getnucharmresults(){

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship =6e+20; //replace to have multiple years of data taking

 TFile *nuefile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_e_dis_cc_charm.root");
 TFile *numufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_mu_dis_cc_charm.root");
 TFile *nutaufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_tau_dis_cc_charm.root");
 TFile *nuebarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_e_bar_dis_cc_charm.root");
 TFile *numubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_mu_bar_dis_cc_charm.root");
 TFile *nutaubarfile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/plots_2/results_nu_tau_bar_dis_cc_charm.root");


 TH1D *nuespectrum = (TH1D*) nuefile->Get("hspectrum_nu_e_intdis_cc_charm");
 TH1D *numuspectrum = (TH1D*) numufile->Get("hspectrum_nu_mu_intdis_cc_charm");
 TH1D *nutauspectrum = (TH1D*) nutaufile->Get("hspectrum_nu_tau_intdis_cc_charm");

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

 TH1D *nuebarspectrum = (TH1D*) nuebarfile->Get("hspectrum_nu_e_bar_intdis_cc_charm");
 TH1D *numubarspectrum = (TH1D*) numubarfile->Get("hspectrum_nu_mu_bar_intdis_cc_charm");
 TH1D *nutaubarspectrum = (TH1D*) nutaubarfile->Get("hspectrum_nu_tau_bar_intdis_cc_charm");

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
 double normship = 5e+13; //replace to have multiple years of data taking

 TFile *nufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/neutrinos_detector.root");

 TH1D *nuespectrum = (TH1D*) nufile->Get("hnu_e");
 TH1D *numuspectrum = (TH1D*) nufile->Get("hnu_mu");
 TH1D *nutauspectrum = (TH1D*) nufile->Get("hnu_tau");

 nuespectrum->Scale(normship/normsim);
 numuspectrum->Scale(normship/normsim);
 nutauspectrum->Scale(normship/normsim);
 //printing yields to screen
 cout<<"nue: Mean energy<<"<<nuespectrum->GetMean()<<" N:"<<nuespectrum->Integral()<<endl;
 cout<<endl;
 cout<<"numu: Mean energy<<"<<numuspectrum->GetMean()<<" N:"<<numuspectrum->Integral()<<endl;
 cout<<endl;
 cout<<"nutau: Mean energy<<"<<nutauspectrum->GetMean()<<" N:"<<nutauspectrum->Integral()<<endl;
 cout<<endl;

 TH1D *nuebarspectrum = (TH1D*) nufile->Get("hnu_e_bar");
 TH1D *numubarspectrum = (TH1D*) nufile->Get("hnu_mu_bar");
 TH1D *nutaubarspectrum = (TH1D*) nufile->Get("hnu_tau_bar");

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

void getnuproducedresults(){
 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 6e+20;//replace to have multiple years of data taking

 TFile *nufile = TFile::Open("/home/utente/Simulations/nuyield_shipecn3/25m/pythia8_Geant4_1.0_withCharm_nu.root");

 TH1D *nuespectrum = (TH1D*) nufile->Get("1012");
 TH1D *numuspectrum = (TH1D*) nufile->Get("1014");
 TH1D *nutauspectrum = (TH1D*) nufile->Get("1016");

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

 TH1D *nuebarspectrum = (TH1D*) nufile->Get("2012");
 TH1D *numubarspectrum = (TH1D*) nufile->Get("2014");
 TH1D *nutaubarspectrum = (TH1D*) nufile->Get("2016");

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

 //combining nu and antinu together
 nuespectrum->Add(nuebarspectrum);
 numuspectrum->Add(numubarspectrum);
 nutauspectrum->Add(nutaubarspectrum);

 nuespectrum->SetTitle("#nu_e + anti-#nu_e");
 numuspectrum->SetTitle("#nu_#mu + anti-#nu_#mu");
 nutauspectrum->SetTitle("#nu_#tau + anti-#nu_#tau");

 TCanvas *cnu = new TCanvas();
 numuspectrum->SetLineColor(kRed);
 numuspectrum->Draw();
 nuespectrum->SetLineColor(kBlue);
 nuespectrum->Draw("SAMES");
 nutauspectrum->SetLineColor(kGreen);
 nutauspectrum->Draw("SAMES");

 //numuspectrum->GetYaxis()->SetRangeUser(1e+3,3e+10);
 cnu->BuildLegend();

}