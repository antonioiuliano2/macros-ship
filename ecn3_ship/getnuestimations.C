//just draw histograms together, with colors as in SHiP CDS plots, get them together
//TString prepath("/home/utente/Simulations/nuyield_shipecn3/advsnd/");
//TString prepath("/home/utente/Simulations/nuyield_shipecn3/advsnd_downstream/");
TString prepath("/home/utente/Simulations/nuyield_shipecn3/2025_08_28_nuyield_SND_EmuTargetSiliconTarget/");
void getnuresults(){

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 4e+19; //replace to have multiple years of data taking

 TFile *nuefile = TFile::Open(prepath+TString("plots_2/results_nu_e_dis_cc.root").Data());
 TFile *numufile = TFile::Open(prepath+TString("plots_2/results_nu_mu_dis_cc.root").Data());
 TFile *nutaufile = TFile::Open(prepath+TString("plots_2/results_nu_tau_dis_cc.root").Data());
 TFile *nuebarfile = TFile::Open(prepath+TString("plots_2/results_nu_e_bar_dis_cc.root").Data());
 TFile *numubarfile = TFile::Open(prepath+TString("plots_2/results_nu_mu_bar_dis_cc.root").Data());
 TFile *nutaubarfile = TFile::Open(prepath+TString("plots_2/results_nu_tau_bar_dis_cc.root").Data());


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

 nuespectrum->SetTitle("#nu_{e} + anti#nu_{e}");
 numuspectrum->SetTitle("#nu_{#mu} + anti#nu_{#mu}");
 nutauspectrum->SetTitle("#nu_{#tau} + anti#nu_{#tau}");
 cnu->BuildLegend();


}

void getnucharmresults(){

 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship =4e+19; //replace to have multiple years of data taking

 TFile *nuefile = TFile::Open(prepath+TString("plots_2/results_nu_e_dis_cc_charm.root").Data());
 TFile *numufile = TFile::Open(prepath+TString("plots_2/results_nu_mu_dis_cc_charm.root").Data());
 TFile *nutaufile = TFile::Open(prepath+TString("plots_2/results_nu_tau_dis_cc_charm.root").Data());
 TFile *nuebarfile = TFile::Open(prepath+TString("plots_2/results_nu_e_bar_dis_cc_charm.root").Data());
 TFile *numubarfile = TFile::Open(prepath+TString("plots_2/results_nu_mu_bar_dis_cc_charm.root").Data());
 TFile *nutaubarfile = TFile::Open(prepath+TString("plots_2/results_nu_tau_bar_dis_cc_charm.root").Data());


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

 nuespectrum->SetTitle("nu_e + anti-nu_e");
 numuspectrum->SetTitle("nu_mu + anti-nu_mu");
 nutauspectrum->SetTitle("nu_tau + anti-nu_tau");
 cnu->BuildLegend();


}

void getnuincomingresults(){
 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 4e+19; //replace to have multiple years of data taking

 TFile *nufile = TFile::Open(prepath+TString("neutrinos_detector.root").Data());

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

 nuespectrum->SetTitle("nu_e + anti-nu_e");
 numuspectrum->SetTitle("nu_mu + anti-nu_mu");
 nutauspectrum->SetTitle("nu_tau + anti-nu_tau");
 cnu->BuildLegend();

}

void getnuproducedresults(){
 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 4e+19;//replace to have multiple years of data taking

 TFile *nufile = TFile::Open(prepath+TString("pythia8_Geant4_1.0_withCharm_nu.root").Data());

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

 nuespectrum->SetTitle("#nu_{e} + anti#nu_{e}");
 numuspectrum->SetTitle("#nu_{#mu} + anti#nu_{#mu}");
 nutauspectrum->SetTitle("#nu_{#tau} + anti#nu_{#tau}");

 TCanvas *cnu = new TCanvas();
 numuspectrum->GetXaxis()->SetTitle("GeV/c");
 numuspectrum->SetLineColor(kRed);
 numuspectrum->Draw("histo");
 nuespectrum->GetXaxis()->SetTitle("GeV/c");
 nuespectrum->SetLineColor(kBlue);
 nuespectrum->Draw("histo&&SAMES");
 nutauspectrum->GetXaxis()->SetTitle("GeV/c");
 nutauspectrum->SetLineColor(kGreen);
 nutauspectrum->Draw("histo&&SAMES");

 //numuspectrum->GetYaxis()->SetRangeUser(1e+3,3e+10);
 cnu->BuildLegend();

}

//************************************SAME AS BEFORE, but comparing different production modes, ***********************//
/*
 From Thomas README:
 pythia8_Geant4_1.0_c_nu.root (mbias without charm)
 pythia8_Geant4_charm_nu_1.0.root (charm from cascade)
 pythia8_Geant4_1.0_withCharm_nu.root  (sum of above)

*/

void compareparent_nuproduced(){
 //charm from cascade
 TFile *nufile_charm = TFile::Open(prepath+TString("pythia8_Geant4_charm_nu_1.0.root").Data());

 TH1D *nuespectrum_charm = (TH1D*) nufile_charm->Get("12");
 TH1D *numuspectrum_charm = (TH1D*) nufile_charm->Get("14");
 TH1D *nutauspectrum_charm = (TH1D*) nufile_charm->Get("16");
 
 TH1D *nuebarspectrum_charm = (TH1D*) nufile_charm->Get("-12");
 TH1D *numubarspectrum_charm = (TH1D*) nufile_charm->Get("-14");
 TH1D *nutaubarspectrum_charm = (TH1D*) nufile_charm->Get("-16");

 //combining nu and antinu together
 nuespectrum_charm->Add(nuebarspectrum_charm);
 numuspectrum_charm->Add(numubarspectrum_charm);
 nutauspectrum_charm->Add(nutaubarspectrum_charm);

 nuespectrum_charm->SetName("nu_e_charm");
 numuspectrum_charm->SetName("nu_mu_charm");
 nutauspectrum_charm->SetName("nu_tau_charm");

 nuespectrum_charm->SetTitle("nu_e + anti-nu_e charm from cascade");
 numuspectrum_charm->SetTitle("nu_mu + anti-nu_mu charm from cascade");
 nutauspectrum_charm->SetTitle("nu_tau + anti-nu_tau charm from cascade");

 //(mbias without charm)

 TFile *nufile_c = TFile::Open(prepath+TString("pythia8_Geant4_1.0_c_nu.root").Data());

 TH1D *nuespectrum_c = (TH1D*) nufile_c->Get("12");
 TH1D *numuspectrum_c = (TH1D*) nufile_c->Get("14");
 TH1D *nutauspectrum_c = (TH1D*) nufile_c->Get("16");
 
 TH1D *nuebarspectrum_c = (TH1D*) nufile_c->Get("-12");
 TH1D *numubarspectrum_c = (TH1D*) nufile_c->Get("-14");
 TH1D *nutaubarspectrum_c = (TH1D*) nufile_c->Get("-16");

 //combining nu and antinu together
 nuespectrum_c->Add(nuebarspectrum_c);
 numuspectrum_c->Add(numubarspectrum_c);
 nutauspectrum_c->Add(nutaubarspectrum_c);

 const double total_nue = nuespectrum_c->Integral() + nuespectrum_charm->Integral();
 const double total_numu = numuspectrum_c->Integral() + numuspectrum_charm->Integral();
 const double total_nutau = nutauspectrum_c->Integral() + nutauspectrum_charm->Integral();

 //scaling to the total number of neutrinos, both from mbias and charm
 nuespectrum_charm->Scale(1./total_nue);
 numuspectrum_charm->Scale(1./total_numu);
 nutauspectrum_charm->Scale(1./total_nutau);

 nuespectrum_c->Scale(1./total_nue);
 numuspectrum_c->Scale(1./total_numu);
 nutauspectrum_c->Scale(1./total_nutau);

 nuespectrum_c->SetName("nu_e_mbias");
 numuspectrum_c->SetName("nu_mu_mbias");
 nutauspectrum_c->SetName("nu_tau_mbias");
 nuespectrum_c->SetTitle("nu_e + anti-nu_e mbias without charm");
 numuspectrum_c->SetTitle("nu_mu + anti-nu_mu mbias without charm");
 nutauspectrum_c->SetTitle("nu_tau + anti-nu_tau mbias without charm");

 //plotting distributions
 TCanvas *cnu_e = new TCanvas();
 nuespectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nuespectrum_charm->SetLineColor(kRed);
 nuespectrum_c->SetLineColor(kBlue);
 nuespectrum_c->Draw("histo");
 nuespectrum_charm->Draw("SAMES&&histo");
 cnu_e->SetLogy();
 cnu_e->BuildLegend();

 TCanvas *cnu_mu = new TCanvas();
 numuspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 numuspectrum_charm->SetLineColor(kRed);
 numuspectrum_c->SetLineColor(kBlue);
 numuspectrum_c->Draw("histo");
 numuspectrum_charm->Draw("SAMES&&histo");
 cnu_mu->SetLogy();
 cnu_mu->BuildLegend();

 TCanvas *cnu_tau = new TCanvas();
 nutauspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nutauspectrum_charm->SetLineColor(kRed);
 nutauspectrum_charm->Draw("histo");
 nutauspectrum_c->SetLineColor(kBlue);
 nutauspectrum_c->Draw("SAMES&&histo");
 cnu_tau->SetLogy();
 cnu_tau->BuildLegend();

 gStyle->SetOptStat("nei");

 cout<<"Nue Fractions "<<nuespectrum_c->Integral()<<" "<<nuespectrum_charm->Integral()<<endl;
 cout<<"Numu Fractions "<<numuspectrum_c->Integral()<<" "<<numuspectrum_charm->Integral()<<endl;
 cout<<"Nutau Fractions "<<nutauspectrum_c->Integral()<<" "<<nutauspectrum_charm->Integral()<<endl;

}


void compareparent_nudetector(){
 //charm from cascade
 TFile *nufile_charm = TFile::Open(prepath+TString("neutrinos_detector_charm.root").Data());

 TH1D *nuespectrum_charm = (TH1D*) nufile_charm->Get("hnu_e");
 TH1D *numuspectrum_charm = (TH1D*) nufile_charm->Get("hnu_mu");
 TH1D *nutauspectrum_charm = (TH1D*) nufile_charm->Get("hnu_tau");
 
 TH1D *nuebarspectrum_charm = (TH1D*) nufile_charm->Get("hnu_e_bar");
 TH1D *numubarspectrum_charm = (TH1D*) nufile_charm->Get("hnu_mu_bar");
 TH1D *nutaubarspectrum_charm = (TH1D*) nufile_charm->Get("hnu_tau_bar");

 //combining nu and antinu together
 nuespectrum_charm->Add(nuebarspectrum_charm);
 numuspectrum_charm->Add(numubarspectrum_charm);
 nutauspectrum_charm->Add(nutaubarspectrum_charm);

 nuespectrum_charm->SetName("nu_e_charm");
 numuspectrum_charm->SetName("nu_mu_charm");
 nutauspectrum_charm->SetName("nu_tau_charm");

 nuespectrum_charm->SetTitle("nu_e + anti-nu_e charm from cascade");
 numuspectrum_charm->SetTitle("nu_mu + anti-nu_mu charm from cascade");
 nutauspectrum_charm->SetTitle("nu_tau + anti-nu_tau charm from cascade");

 //(mbias without charm)

 TFile *nufile_c = TFile::Open(prepath+TString("neutrinos_detector_nocharm.root").Data());

 TH1D *nuespectrum_c = (TH1D*) nufile_c->Get("hnu_e");
 TH1D *numuspectrum_c = (TH1D*) nufile_c->Get("hnu_mu");
 TH1D *nutauspectrum_c = (TH1D*) nufile_c->Get("hnu_tau");
 
 TH1D *nuebarspectrum_c = (TH1D*) nufile_c->Get("hnu_e_bar");
 TH1D *numubarspectrum_c = (TH1D*) nufile_c->Get("hnu_mu_bar");
 TH1D *nutaubarspectrum_c = (TH1D*) nufile_c->Get("hnu_tau_bar");

 //combining nu and antinu together
 nuespectrum_c->Add(nuebarspectrum_c);
 numuspectrum_c->Add(numubarspectrum_c);
 nutauspectrum_c->Add(nutaubarspectrum_c);

 const double total_nue = nuespectrum_c->Integral() + nuespectrum_charm->Integral();
 const double total_numu = numuspectrum_c->Integral() + numuspectrum_charm->Integral();
 const double total_nutau = nutauspectrum_c->Integral() + nutauspectrum_charm->Integral();

 //scaling to the total number of neutrinos, both from mbias and charm
 nuespectrum_charm->Scale(1./total_nue);
 numuspectrum_charm->Scale(1./total_numu);
 nutauspectrum_charm->Scale(1./total_nutau);

 nuespectrum_c->Scale(1./total_nue);
 numuspectrum_c->Scale(1./total_numu);
 nutauspectrum_c->Scale(1./total_nutau);

 nuespectrum_c->SetName("nu_e_mbias");
 numuspectrum_c->SetName("nu_mu_mbias");
 nutauspectrum_c->SetName("nu_tau_mbias");
 nuespectrum_c->SetTitle("nu_e + anti-nu_e mbias without charm");
 numuspectrum_c->SetTitle("nu_mu + anti-nu_mu mbias without charm");
 nutauspectrum_c->SetTitle("nu_tau + anti-nu_tau mbias without charm");

 //plotting distributions
 TCanvas *cnu_e = new TCanvas();
 nuespectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nuespectrum_charm->SetLineColor(kRed);
 nuespectrum_c->SetLineColor(kBlue);
 nuespectrum_c->Draw("histo");
 nuespectrum_charm->Draw("SAMES&&histo");
 cnu_e->SetLogy();
 cnu_e->BuildLegend();

 TCanvas *cnu_mu = new TCanvas();
 numuspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 numuspectrum_charm->SetLineColor(kRed);
 numuspectrum_c->SetLineColor(kBlue);
 numuspectrum_c->Draw("histo");
 numuspectrum_charm->Draw("SAMES&&histo");
 cnu_mu->SetLogy();
 cnu_mu->BuildLegend();

 TCanvas *cnu_tau = new TCanvas();
 nutauspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nutauspectrum_charm->SetLineColor(kRed);
 nutauspectrum_charm->Draw("histo");
 nutauspectrum_c->SetLineColor(kBlue);
 nutauspectrum_c->Draw("SAMES&&histo");
 cnu_tau->SetLogy();
 cnu_tau->BuildLegend();

 gStyle->SetOptStat("nei");

 cout<<"Nue Fractions "<<nuespectrum_c->Integral()<<" "<<nuespectrum_charm->Integral()<<endl;
 cout<<"Numu Fractions "<<numuspectrum_c->Integral()<<" "<<numuspectrum_charm->Integral()<<endl;
 cout<<"Nutau Fractions "<<nutauspectrum_c->Integral()<<" "<<nutauspectrum_charm->Integral()<<endl;

}

void compareparent_nuinteracting(){
 //charm from cascade
 TFile *nuefile_charm = TFile::Open(prepath+TString("plots_1/results_nu_e_dis_cc.root").Data());
 TH1D *nuespectrum_charm = (TH1D*) nuefile_charm->Get("hspectrum_nu_e_intdis_cc");

 TFile *numufile_charm = TFile::Open(prepath+TString("plots_1/results_nu_mu_dis_cc.root").Data());
 TH1D *numuspectrum_charm = (TH1D*) numufile_charm->Get("hspectrum_nu_mu_intdis_cc");

 TFile *nutaufile_charm = TFile::Open(prepath+TString("plots_1/results_nu_tau_dis_cc.root").Data());
 TH1D *nutauspectrum_charm = (TH1D*) nutaufile_charm->Get("hspectrum_nu_tau_intdis_cc");

 TFile *nue_bar_file_charm = TFile::Open(prepath+TString("plots_1/results_nu_e_bar_dis_cc.root").Data());
 TH1D *nuebarspectrum_charm = (TH1D*) nue_bar_file_charm->Get("hspectrum_nu_e_bar_intdis_cc");

 TFile *numu_bar_file_charm = TFile::Open(prepath+TString("plots_1/results_nu_mu_bar_dis_cc.root").Data());
 TH1D *numubarspectrum_charm = (TH1D*) numu_bar_file_charm->Get("hspectrum_nu_mu_bar_intdis_cc");

 TFile *nutau_bar_file_charm = TFile::Open(prepath+TString("plots_1/results_nu_tau_bar_dis_cc.root").Data());
 TH1D *nutaubarspectrum_charm = (TH1D*) nutau_bar_file_charm->Get("hspectrum_nu_tau_bar_intdis_cc");

 //combining nu and antinu together
 nuespectrum_charm->Add(nuebarspectrum_charm);
 numuspectrum_charm->Add(numubarspectrum_charm);
 nutauspectrum_charm->Add(nutaubarspectrum_charm);

 nuespectrum_charm->SetName("nu_e_charm");
 numuspectrum_charm->SetName("nu_mu_charm");
 nutauspectrum_charm->SetName("nu_tau_charm");

 nuespectrum_charm->SetTitle("nu_e + anti-nu_e charm from cascade");
 numuspectrum_charm->SetTitle("nu_mu + anti-nu_mu charm from cascade");
 nutauspectrum_charm->SetTitle("nu_tau + anti-nu_tau charm from cascade");

 //(mbias without charm)

 TFile *nuefile_c = TFile::Open(prepath+TString("plots_0/results_nu_e_dis_cc.root").Data());
 TH1D *nuespectrum_c = (TH1D*) nuefile_c->Get("hspectrum_nu_e_intdis_cc");

 TFile *numufile_c = TFile::Open(prepath+TString("plots_0/results_nu_mu_dis_cc.root").Data());
 TH1D *numuspectrum_c = (TH1D*) numufile_c->Get("hspectrum_nu_mu_intdis_cc");

 TFile *nutaufile_c = TFile::Open(prepath+TString("plots_0/results_nu_tau_dis_cc.root").Data());
 TH1D *nutauspectrum_c = (TH1D*) nutaufile_c->Get("hspectrum_nu_tau_intdis_cc");

 TFile *nue_bar_file_c = TFile::Open(prepath+TString("plots_0/results_nu_e_bar_dis_cc.root").Data());
 TH1D *nuebarspectrum_c = (TH1D*) nue_bar_file_c->Get("hspectrum_nu_e_bar_intdis_cc");

 TFile *numu_bar_file_c = TFile::Open(prepath+TString("plots_0/results_nu_mu_bar_dis_cc.root").Data());
 TH1D *numubarspectrum_c = (TH1D*) numu_bar_file_c->Get("hspectrum_nu_mu_bar_intdis_cc");

 TFile *nutau_bar_file_c = TFile::Open(prepath+TString("plots_0/results_nu_tau_bar_dis_cc.root").Data());
 TH1D *nutaubarspectrum_c = (TH1D*) nutau_bar_file_c->Get("hspectrum_nu_tau_bar_intdis_cc");

 //combining nu and antinu together
 nuespectrum_c->Add(nuebarspectrum_c);
 numuspectrum_c->Add(numubarspectrum_c);
 nutauspectrum_c->Add(nutaubarspectrum_c);

 const double total_nue = nuespectrum_c->Integral() + nuespectrum_charm->Integral();
 const double total_numu = numuspectrum_c->Integral() + numuspectrum_charm->Integral();
 const double total_nutau = nutauspectrum_c->Integral() + nutauspectrum_charm->Integral();

 cout<<"Nue Fractions "<<nuespectrum_c->Integral()<<" "<<nuespectrum_charm->Integral()<<endl;
 cout<<"Numu Fractions "<<numuspectrum_c->Integral()<<" "<<numuspectrum_charm->Integral()<<endl;
 cout<<"Nutau Fractions "<<nutauspectrum_c->Integral()<<" "<<nutauspectrum_charm->Integral()<<endl;

 //scaling to the total number of neutrinos, both from mbias and charm
 nuespectrum_charm->Scale(1./total_nue);
 numuspectrum_charm->Scale(1./total_numu);
 nutauspectrum_charm->Scale(1./total_nutau);

 nuespectrum_c->Scale(1./total_nue);
 numuspectrum_c->Scale(1./total_numu);
 nutauspectrum_c->Scale(1./total_nutau);

 nuespectrum_c->SetName("nu_e_mbias");
 numuspectrum_c->SetName("nu_mu_mbias");
 nutauspectrum_c->SetName("nu_tau_mbias");
 nuespectrum_c->SetTitle("nu_e + anti-nu_e mbias without charm");
 numuspectrum_c->SetTitle("nu_mu + anti-nu_mu mbias without charm");
 nutauspectrum_c->SetTitle("nu_tau + anti-nu_tau mbias without charm");

 //plotting distributions
 TCanvas *cnu_e = new TCanvas();
 nuespectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nuespectrum_charm->SetLineColor(kRed);
 nuespectrum_c->SetLineColor(kBlue);
 nuespectrum_charm->Draw("histo");
 nuespectrum_c->Draw("SAMES&&histo");
 cnu_e->SetLogy();
 cnu_e->BuildLegend();

 TCanvas *cnu_mu = new TCanvas();
 numuspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 numuspectrum_charm->SetLineColor(kRed);
 numuspectrum_c->SetLineColor(kBlue);
 numuspectrum_c->Draw("histo");
 numuspectrum_charm->Draw("SAMES&&histo");
 cnu_mu->SetLogy();
 cnu_mu->BuildLegend();

 TCanvas *cnu_tau = new TCanvas();
 nutauspectrum_charm->GetXaxis()->SetTitle("p[GeV/c]");
 nutauspectrum_charm->SetLineColor(kRed);
 nutauspectrum_charm->Draw("histo");
 nutauspectrum_c->SetLineColor(kBlue);
 nutauspectrum_c->Draw("SAMES&&histo");
 cnu_tau->SetLogy();
 cnu_tau->BuildLegend();

 gStyle->SetOptStat("nei");

 cout<<"Nue Fractions "<<nuespectrum_c->Integral()<<" "<<nuespectrum_charm->Integral()<<endl;
 cout<<"Numu Fractions "<<numuspectrum_c->Integral()<<" "<<numuspectrum_charm->Integral()<<endl;
 cout<<"Nutau Fractions "<<nutauspectrum_c->Integral()<<" "<<nutauspectrum_charm->Integral()<<endl;

}