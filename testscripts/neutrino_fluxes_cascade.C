#include <map>
//projecting neutrino spectra at detector z, now separately for primary and cascade charm production (15 Giugno 2020)
void neutrino_fluxes_cascade(){ //projecting neutrino fluxes to the target
 TFile * f = TFile::Open("Decay-Cascade10M-parp16-MSTP82-1-MSEL4-ntuple.root");
 TTree *tree = (TTree*) f->Get("Decay");

 Double_t start[3]; //start and end of neutrino positions at z of neutrino target
 const Int_t nneutrinos = tree->GetEntries();
 //getting the branches

 Float_t id, px, py, pz,w,x,y;
 Float_t Cascadek;

 tree->SetBranchAddress("id",&id);
 tree->SetBranchAddress("px",&px);
 tree->SetBranchAddress("py",&py);
 tree->SetBranchAddress("pz",&pz);
 tree->SetBranchAddress("weight",&w);
 tree->SetBranchAddress("Cascadek",&Cascadek);
 
 Float_t deltaz = 3969.; //distance between center of proton target and start of neutrino target
 //Float_t deltaz = 5500.;
 Float_t tanx,tany; //angles and positions after projections

 Float_t nnudet = 0., nallneutrinos = 0.;
 map<Int_t, Float_t> nall = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos produced, mapped per pdg
 map<Int_t, Float_t> ndet = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos arrived at det

 //now need to fill differently primary and cascade
 map<Int_t, TH1D*> hspectrumdet_primary;
 map<Int_t, TH1D*> hspectrumdet_cascade;

 TFile *outfile = new TFile("neutrinos_detector_decays10M.root","RECREATE"); //spectra arrived in detector are saved here
 //nutau
 hspectrumdet_primary[16] = new TH1D("hnu_tau_primary","Spectrum tau neutrinos arrived at detector",400,0,400);
 hspectrumdet_primary[-16] = new TH1D("hnu_tau_bar_primary","Spectrum tau neutrinos arrived at detector",400,0,400);

 hspectrumdet_cascade[16] = new TH1D("hnu_tau_cascade","Spectrum tau neutrinos arrived at detector",400,0,400);
 hspectrumdet_cascade[-16] = new TH1D("hnu_tau_bar_cascade","Spectrum tau neutrinos arrived at detector",400,0,400);

 //numu
 hspectrumdet_primary[14] = new TH1D("hnu_mu_primary","Spectrum muon neutrinos arrived at detector",400,0,400);
 hspectrumdet_primary[-14] = new TH1D("hnu_mu_bar_primary","Spectrum muon antineutrinos arrived at detector",400,0,400);

 hspectrumdet_cascade[14] = new TH1D("hnu_mu_cascade","Spectrum muon neutrinos arrived at detector",400,0,400);
 hspectrumdet_cascade[-14] = new TH1D("hnu_mu_bar_cascade","Spectrum muon antineutrinos arrived at detector",400,0,400);

 //nue

 hspectrumdet_primary[12] = new TH1D("hnu_e_primary","Spectrum electron neutrinos arrived at detector",400,0,400);
 hspectrumdet_primary[-12] = new TH1D("hnu_e_bar_primary","Spectrum electron antineutrinos arrived at detector",400,0,400);

 hspectrumdet_cascade[12] = new TH1D("hnu_e_cascade","Spectrum electron neutrinos arrived at detector",400,0,400);
 hspectrumdet_cascade[-12] = new TH1D("hnu_e_bar_cascade","Spectrum electron antineutrinos arrived at detector",400,0,400);
 
 Double_t targetdx = 40., targetdy = 40.;
 //Double_t targetdx = 55, targetdy = 55;
 cout<<"N NEUTRINOS AND MUONS: "<<nneutrinos<<endl;
 for (int i = 0; i < nneutrinos; i++){

 tree->GetEntry(i);
 if (TMath::Abs(id)==13) continue;
 if (i%100000 == 0) cout<<i<<endl;
 tanx = px/pz;
 tany = py/pz;

 Double_t momentum = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

 start[0] = tanx * deltaz; //projecting produced neutrinos to neutrino detector, aggiungere x e y non cambia significativamente il risultato
 start[1] = tany * deltaz;
 start[2] = -3259;

 nallneutrinos += w;

 nall[id] +=w;

 if(TMath::Abs(start[0]) < targetdx && TMath::Abs(start[1]) < targetdy){ //checking how many neutrinos are inside the detector
  nnudet += w;
  ndet[id] += w;
  if(Cascadek==1) hspectrumdet_primary[id]->Fill(momentum,w);
  else hspectrumdet_cascade[id]->Fill(momentum,w);
  }
 }
 
// cout<<hxy_all->Integral()<<endl; //totale entries
 cout<<nallneutrinos<<endl;
 cout<<nnudet<<endl; //dovrebbe escludere underflow and overflow
 cout<<"Ratio arriving at det: "<<nnudet/nallneutrinos<<endl;

 cout<<"Electron neutrinos"<<endl;
 cout<<nall[12]<<" "<<nall[-12]<<endl;
 cout<<ndet[12]<<" "<<ndet[-12]<<endl; //dovrebbe escludere underflow and overflow
 cout<<"Ratio arriving at det: "<<ndet[12]/nall[12]<<" "<<ndet[-12]/nall[-12]<<endl;;

 cout<<"Muon neutrinos"<<endl;
 cout<<nall[14]<<" "<<nall[-14]<<endl;
 cout<<ndet[14]<<" "<<ndet[-14]<<" "<<endl; //dovrebbe escludere underflow and overflow
 cout<<"Ratio arriving at det: "<<ndet[14]/nall[14]<<" "<<ndet[-14]/nall[-14]<<endl;

 cout<<"Tau neutrinos and anti neutrinos"<<endl;
 cout<<nall[16]<<" "<<nall[-16]<<endl;
 cout<<ndet[16]<<" "<<ndet[-16]<<endl; //dovrebbe escludere underflow and overflow
 cout<<"Ratio arriving at det: "<<ndet[16]/nall[16]<<" "<<ndet[-16]/nall[-16]<<endl;

 outfile->Write();
 outfile->Close();
}

void drawhistos(){
 //draw obtained spectra
 TFile *producedspectra = TFile::Open("Decay-Cascade10M-parp16-MSTP82-1-MSEL4-ntuple.root");
 TFile *projectedspectra = TFile::Open("neutrinos_detector_decays10M.root");

 float fDs = 0.077 / 0.106; //correct for too many Ds produced by Pythia6 

 TH1D* hnutau_primary = (TH1D*) producedspectra->Get("1016");
 TH1D* hnutau_cascade = (TH1D*) producedspectra->Get("1316");
 hnutau_cascade->SetLineColor(kRed);

 TH1D* hnutau_projprimary = (TH1D*) projectedspectra->Get("hnu_tau_primary");
 TH1D* hnutau_projcascade = (TH1D*) projectedspectra->Get("hnu_tau_cascade");
 hnutau_projcascade->SetLineColor(kRed);

 TH1D* hnutaubar_primary = (TH1D*) producedspectra->Get("2016");
 TH1D* hnutaubar_cascade = (TH1D*) producedspectra->Get("2316");
 hnutau_cascade->SetLineColor(kRed);

 TH1D* hnutaubar_projprimary = (TH1D*) projectedspectra->Get("hnu_tau_bar_primary");
 TH1D* hnutaubar_projcascade = (TH1D*) projectedspectra->Get("hnu_tau_bar_cascade");
 hnutau_projcascade->SetLineColor(kRed);

 gStyle->SetOptStat(1111);


 hnutau_primary->Scale(fDs);
 hnutau_cascade->Scale(fDs);
 hnutau_projprimary->Scale(fDs);
 hnutau_projcascade->Scale(fDs);

 hnutaubar_primary->Scale(fDs);
 hnutaubar_cascade->Scale(fDs);
 hnutaubar_projprimary->Scale(fDs);
 hnutaubar_projcascade->Scale(fDs);

 //report number neutrinos
 cout<<"In one spill of 5e+15 protons"<<endl;
 cout<<"# tau neutrinos are produced: "<<hnutau_primary->Integral()<<" directly  and "<<hnutau_cascade->Integral()<<" from cascade total: "<< hnutau_primary->Integral() + hnutau_cascade->Integral()<<endl;
 cout<<"# tau neutrinos arrive at detector: "<<hnutau_projprimary->Integral()<<" directly  and "<<hnutau_projcascade->Integral()<<" from cascade total: "<<hnutau_projprimary->Integral()+hnutau_projcascade->Integral()<<endl;
 cout<<"# tau antineutrinos are produced: "<<hnutaubar_primary->Integral()<<" directly  and "<<hnutaubar_cascade->Integral()<<" from cascade total: "<<hnutaubar_primary->Integral()+hnutaubar_cascade->Integral()<<endl;
 cout<<"# tau antineutrinos arrive at detector: "<<hnutaubar_projprimary->Integral()<<" directly  and "<<hnutaubar_projcascade->Integral()<<" from cascade total: "<<hnutaubar_projprimary->Integral()+hnutaubar_projcascade->Integral()<<endl;

 //drawing canvases (summing nutau and antinutau)

 TCanvas *cproduced = new TCanvas();
 hnutau_primary->Add(hnutaubar_primary);
 hnutau_cascade->Add(hnutaubar_cascade);
 hnutau_cascade->Draw();
 hnutau_primary->Draw("SAMES");
 hnutau_primary->SetTitle("Produced nutau and antinutau, primary");
 hnutau_cascade->SetTitle("Produced nutau and antinutau, cascade");
 cproduced->SetLogy();
 cproduced->BuildLegend();

 TCanvas *cprojected = new TCanvas();
 hnutau_projprimary->Add(hnutaubar_projprimary);
 hnutau_projcascade->Add(hnutaubar_projcascade);
 hnutau_projprimary->Draw();
 hnutau_projcascade->Draw("SAMES");
 hnutau_projprimary->SetTitle("At detector nutau and antinutau, primary");
 hnutau_projcascade->SetTitle("At detector nutau and antinutau, cascade");
 cprojected->SetLogy();
 cprojected->BuildLegend();
 
 //2d spectra
 TH2D *hnutau_ppt_primary = (TH2D*) producedspectra->Get("1216");
 TH2D *hnutau_ppt_cascade = (TH2D*) producedspectra->Get("1516");

 TCanvas *cppt = new TCanvas();
 hnutau_ppt_cascade->SetMarkerColor(kRed);
 hnutau_ppt_primary->Draw();
 hnutau_ppt_cascade->Draw("SAMES");

 //interacting
 TFile *interactingspectra_nutau = TFile::Open("plots/results_nu_tau_dis_cc.root");
 TFile *interactingspectra_nutaubar = TFile::Open("plots/results_nu_tau_bar_dis_cc.root");

 TH1D *hspectrum_nutau_primary = (TH1D*) interactingspectra_nutau->Get("hspectrum_nu_tau_intdis_cc");
 TH1D *hspectrum_nutau_cascade = (TH1D*) interactingspectra_nutau->Get("hspectrum_nu_tau_intcascadedis_cc");
 TH1D *hspectrum_nutaubar_primary = (TH1D*) interactingspectra_nutaubar->Get("hspectrum_nu_tau_bar_intdis_cc");
 TH1D *hspectrum_nutaubar_cascade = (TH1D*) interactingspectra_nutaubar->Get("hspectrum_nu_tau_bar_intcascadedis_cc");

 //applying fDs correction
 hspectrum_nutau_primary->Scale(fDs);
 hspectrum_nutaubar_primary->Scale(fDs);
 hspectrum_nutau_cascade->Scale(fDs);
 hspectrum_nutaubar_cascade->Scale(fDs);

 //summing nutau and nutaubar
 hspectrum_nutau_primary->Add(hspectrum_nutaubar_primary);
 hspectrum_nutaubar_primary->Add(hspectrum_nutaubar_cascade);

 TCanvas *cdis = new TCanvas();
 hspectrum_nutau_primary->Draw("histo");
 hspectrum_nutau_cascade->SetLineColor(kRed);
 hspectrum_nutau_cascade->Draw("histo &&SAMES");

 cout<<hspectrum_nutau_primary->Integral(0,400)<<" "<<hspectrum_nutau_cascade->Integral(0,400)<<endl;
}

///////////////STARTING METHOD FOR ESTIMATING NEUTRINO YIELDS/////////////////////////////////////////////////////
Double_t nu_yield_general(const char* nu = "nu_tau", const char* intmode = "dis_cc", const char* charmmode = "");
void nu_yield(){ //passing neutrino types and interaction mode to nu_yield_general for estimation  
  TCanvas *cspectratau = new TCanvas();
  cspectratau->Divide(3,3);
  cout<<"Yields per nutau"<<endl;
  cspectratau->cd(1);
  Double_t nutau_dis_cc_charm =nu_yield_general("nu_tau","dis_cc","_charm");
  cspectratau->cd(2);
  Double_t nutau_bar_dis_cc_charm =nu_yield_general("nu_tau_bar","dis_cc","_charm");
  cspectratau->cd(3);
  Double_t nutau_dis_cc = nu_yield_general("nu_tau","dis_cc");
  cspectratau->cd(4);
  Double_t nutau_bar_dis_cc =nu_yield_general("nu_tau_bar","dis_cc");
  cspectratau->cd(5);
  Double_t nutau_qel_cc =nu_yield_general("nu_tau","qel_cc");
  cspectratau->cd(6);
  Double_t nutau_bar_qel_cc =nu_yield_general("nu_tau_bar","qel_cc");
  cspectratau->cd(7);  
  Double_t nutau_res_cc =nu_yield_general("nu_tau","res_cc");
  cspectratau->cd(8);
  Double_t nutau_bar_res_cc =nu_yield_general("nu_tau_bar","res_cc");
  cspectratau->cd(9);
  Double_t nutau_el =nu_yield_general("nu_tau","ve_nc");   
  Double_t nutau_bar_el =nu_yield_general("nu_tau_bar","ve_nc");

  cout<<"NUTAU"<<" "<<nutau_dis_cc<<" "<<nutau_dis_cc_charm<<" "<<nutau_res_cc<<" "<<nutau_qel_cc<<" "<<nutau_el<<endl;
  cout<<"ANTINUTAU"<<" "<<nutau_bar_dis_cc<<" "<<nutau_bar_dis_cc_charm<<" "<<nutau_bar_res_cc<<" "<<nutau_bar_qel_cc<<" "<<nutau_bar_el<<endl;
  
}
//general layout with FORM to estimate number of neutrino interactions 
Double_t nu_yield_general(const char* nu = "nu_mu", const char* intmode = "dis_cc", const char* charmmode = ""){     
  TFile *xsec = TFile::Open("Nu_xsec_full.root"); //normal splines are cut at 350 GeV
  TFile *flux = TFile::Open("neutrinos_detector_decays10M.root");

  //const char* nu = "nu_mu"; //neutrino type
  //const char* intmode = "dis_cc"; //interaction mode;

  //  const char* charmmode = "_charm"; //write _charm for CHARMCCDIS
  
  TH1D *hfluxnu = (TH1D*) flux->Get(Form("h%s_primary",nu));
  TH1D *hfluxnucascade = (TH1D*) flux->Get(Form("h%s_cascade",nu));

  bool nullp = false, nulln = false, nullother = false; //for DIS both are possible, but not for other interactions
  TGraph *hxsec_p = (TGraph*) xsec->Get(Form("%s_Pb208/%s_p%s",nu,intmode,charmmode));   
  TGraph *hxsec_n = (TGraph*) xsec->Get(Form("%s_Pb208/%s_n%s",nu,intmode,charmmode));
  TGraph *hxsec_other = (TGraph*) xsec->Get(Form("%s_Pb208/%s%s",nu,intmode,charmmode));
  if (hxsec_p == NULL) nullp = true;
  if (hxsec_n == NULL) nulln = true;
  if (hxsec_other == NULL) nullother = true;
  
  const Int_t A = 208; //initial approximation, only pb208 considered
  const Int_t Z = 82;

  Float_t Ninteracting = 0.;
  Float_t Ninteractingcascade = 0.;
  Double_t mass = 8183*1e+3; //mass in grams
  //Double_t surface = 90.3 * 74.9; //surface in square centimetres (rectangular configuration)
  Double_t surface = 80. * 80.; //surface in square centimetres (squared configuration)
  //Double_t surface = 1.4e+4;
  
  Double_t avogadro = 6.022e+23;
  Double_t NT = mass/A * avogadro;  
  Double_t x, y ,ysec, ycascade;

  TFile *outputfile = new TFile(Form("plots/results_%s_%s%s.root",nu,intmode,charmmode),"RECREATE");
  TGraph *hxsec_total = new TGraph(); //summing protons and neutrons when both are present
  TH1D *hspectrum_int = new TH1D(Form("hspectrum_%s_int%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), 400, 0, 400); 
  TH1D *hspectrum_intcascade = new TH1D(Form("hspectrum_%s_intcascade%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), 400, 0, 400);
  TH1D *htest = new TH1D(Form("htest_%s_int%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), 400, 0, 400); 

  for (int n = 0; n < hfluxnu->GetNbinsX(); n++){
  
  y = hfluxnu->GetBinContent(n+1);
  ycascade = hfluxnucascade->GetBinContent(n+1);
  x = hfluxnu->GetBinLowEdge(n+1) + hfluxnu->GetBinWidth(n+1)/2.;
  
  if (nullp && nulln && nullother) cout<<Form("ERROR: no cross sections are present for neutrino %s in interaction %s%s",nu,intmode,charmmode)<<endl;
  if (nullp && nulln) ysec = hxsec_other->Eval(x);
  else if (nullp) ysec = hxsec_n->Eval(x);
  else if (nulln) ysec = hxsec_p->Eval(x);
  else ysec = (hxsec_p->Eval(x) + hxsec_n->Eval(x)); //summing proton and neutron cross section

  hxsec_total->SetPoint(n+1,x,ysec);
 
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnu->GetBinWidth(n+1);  
  Ninteractingcascade = ysec *1e-38 * NT * ycascade/surface * hfluxnucascade->GetBinWidth(n+1);  
  
  htest->SetBinContent(n+1,hfluxnu->GetBinWidth(n+1)*y);
  hspectrum_int->SetBinContent(n+1,Ninteracting);

  hspectrum_intcascade->SetBinContent(n+1,Ninteractingcascade);

  }
  Double_t counting = hfluxnu->Integral();

  Double_t yield =  hspectrum_int->Integral(0,400);
  Double_t yieldcascade =  hspectrum_intcascade->Integral(0,400);
  //cout<<Form("numero neutrini %s per spill interagenti in mode %s%s",nu,intmode,charmmode)<<" "<<hspectrum_int->Integral()<<endl;
  hspectrum_int->Draw();
  hspectrum_int->GetXaxis()->SetTitle("GeV/c");
  hspectrum_intcascade->GetXaxis()->SetTitle("GeV/c");
  hspectrum_intcascade->SetLineColor(kRed);
  hspectrum_intcascade->Draw("SAME");
  // hxsec_total->Draw();
  hxsec_total->SetTitle(Form("cross section from GENIE for %s in interaction mode %s%s",nu,intmode,charmmode));
  outputfile->Write();
  hxsec_total->Write("xsec");
  outputfile->Close();
  //cout<<"TEST "<<htest->GetEntries()<<endl;
  //Double_t countingtest = htest->Integral();
  //cout<<"Test conteggio "<<counting<<"integrale manuale "<<countingtest<<endl;
  return yield+yieldcascade;
}






//leaving old version for backup, if I discover the general on does not work for whatever reason














void nu_yield_old(){
 
 TFile *xsec = TFile::Open("Nu_xsec.root");

 //TFile *flux = TFile::Open("input-flux.root");
 //recovering neutrino fluxes

 TFile *flux = TFile::Open("neutrinos_detector.root");
 TH1D *hfluxnutau = (TH1D*) flux->Get("hnu_tau");
 TH1D *hfluxnutaubar = (TH1D*) flux->Get("hnu_tau_bar");
 TH1D *hfluxnumu = (TH1D*) flux->Get("hnu_mu");
 TH1D *hfluxnumubar = (TH1D*) flux->Get("hnu_mu_bar");
 TH1D *hfluxnue = (TH1D*) flux->Get("hnu_e");
 TH1D *hfluxnuebar = (TH1D*) flux->Get("hnu_e_bar");
 
 const Int_t A = 208; //initial approximation, only pb208 considered
 const Int_t Z = 82;

 //recovering cross sections. Just noticed, they are summed for the number of protons and neutrons

 //QE only on neutrons for neutrinos, on protons for antineutrinos)
 TGraph *hxsec_nutau_qel_n = (TGraph*) xsec->Get("nu_tau_Pb208/qel_cc_n");
 TGraph *hxsec_nutaubar_qel_p = (TGraph*) xsec->Get("nu_tau_bar_Pb208/qel_cc_p");
 //DIS 
 TGraph *hxsec_nutaup = (TGraph*) xsec->Get("nu_tau_Pb208/dis_cc_p");
 TGraph *hxsec_nutaun = (TGraph*) xsec->Get("nu_tau_Pb208/dis_cc_n");
 TGraph *hxsec_nutaubarp = (TGraph*) xsec->Get("nu_tau_bar_Pb208/dis_cc_p");
 TGraph *hxsec_nutaubarn = (TGraph*) xsec->Get("nu_tau_bar_Pb208/dis_cc_n");

 TGraph *hxsec_numu_distotalp = (TGraph*) xsec->Get("nu_mu_Pb208/dis_cc_p");
 TGraph *hxsec_numu_distotaln = (TGraph*) xsec->Get("nu_mu_Pb208/dis_cc_n");
 TGraph *hxsec_numubar_distotalp = (TGraph*) xsec->Get("nu_mu_bar_Pb208/dis_cc_p");
 TGraph *hxsec_numubar_distotaln = (TGraph*) xsec->Get("nu_mu_bar_Pb208/dis_cc_n");
 
 
 //CHARMCCDIS mu 
 TGraph *hxsec_numu_charmp = (TGraph*) xsec->Get("nu_mu_Pb208/dis_cc_p_charm");
 TGraph *hxsec_numu_charmn = (TGraph*) xsec->Get("nu_mu_Pb208/dis_cc_n_charm");
 TGraph *hxsec_numubar_charmp = (TGraph*) xsec->Get("nu_mu_bar_Pb208/dis_cc_p_charm");
 TGraph *hxsec_numubar_charmn = (TGraph*) xsec->Get("nu_mu_bar_Pb208/dis_cc_n_charm");

 TGraph *hxsec_nue_distotalp = (TGraph*) xsec->Get("nu_e_Pb208/dis_cc_p");
 TGraph *hxsec_nue_distotaln = (TGraph*) xsec->Get("nu_e_Pb208/dis_cc_n");
 TGraph *hxsec_nuebar_distotalp = (TGraph*) xsec->Get("nu_e_bar_Pb208/dis_cc_p");
 TGraph *hxsec_nuebar_distotaln = (TGraph*) xsec->Get("nu_e_bar_Pb208/dis_cc_n");
 
 //CHARMCCDIS e 
 TGraph *hxsec_nue_charmp = (TGraph*) xsec->Get("nu_e_Pb208/dis_cc_p_charm");
 TGraph *hxsec_nue_charmn = (TGraph*) xsec->Get("nu_e_Pb208/dis_cc_n_charm");
 TGraph *hxsec_nuebar_charmp = (TGraph*) xsec->Get("nu_e_bar_Pb208/dis_cc_p_charm");
 TGraph *hxsec_nuebar_charmn = (TGraph*) xsec->Get("nu_e_bar_Pb208/dis_cc_n_charm");
 
 //obtain cross sections summing p and n
 TGraph *hxsec_distotal_nutau = new TGraph();
 hxsec_distotal_nutau->SetName("dis_cc_tot_nutau");
 TGraph *hxsec_distotal_nutaubar = new TGraph();
 hxsec_distotal_nutaubar->SetName("dis_cc_tot_nutaubar");

 TGraph *hxsec_distotal_numu = new TGraph();
 hxsec_distotal_numu->SetName("dis_cc_tot_numu");
 TGraph *hxsec_distotal_numubar = new TGraph();
 hxsec_distotal_numubar->SetName("dis_cc_tot_numubar");
 
 TGraph *hxsec_charmtotal_numu = new TGraph();
 hxsec_charmtotal_numu->SetName("charm_tot_numu");
 TGraph *hxsec_charmtotal_numubar = new TGraph();
 hxsec_charmtotal_numubar->SetName("charm_tot_numubar");

 TGraph *hxsec_distotal_nue = new TGraph();
 hxsec_distotal_nue->SetName("dis_cc_tot_nue");
 TGraph *hxsec_distotal_nuebar = new TGraph();
 hxsec_distotal_nuebar->SetName("dis_cc_tot_nuebar");
 
 TGraph *hxsec_charmtotal_nue = new TGraph();
 hxsec_charmtotal_nue->SetName("charm_tot_nue");
 TGraph *hxsec_charmtotal_nuebar = new TGraph();
 hxsec_charmtotal_nuebar->SetName("charm_tot_nuebar");

 TF1 *hxsec = new TF1("operaxsec","(-0.75+0.175301*x+0.0123593*x*x - 0.000185362*x*x*x + 1.7173E-06*x*x*x*x -9.84985E-09*x*x*x*x*x + 3.39435E-11 *x*x*x*x*x*x -6.41045E-14*x*x*x*x*x*x*x +5.07931E-17*x*x*x*x*x*x*x*x)*208",3.6,300);
 hxsec->Draw();

 hfluxnutau->Draw();
 
 //histograms where entries will be saved

 TH1D *hspectrum_nutau = new TH1D("hspectrum_nutau","Spettro neutrini interagenti DIS", 400, 0, 400);
 TH1D *hspectrum_nutaubar = new TH1D("hspectrum_nutaubar","Spettro neutrini interagenti DIS", 400, 0, 400);

 TH1D *hspectrum_nutau_qel = new TH1D("hspectrum_nutau_qel","Spettro neutrini interagenti QEL", 400, 0, 400);
 TH1D *hspectrum_nutaubar_qel = new TH1D("hspectrum_nutaubar_qel","Spettro neutrini interagenti QEL", 400, 0, 400);

 TH1D *hspectrum_nue_charm = new TH1D("hspectrum_nue_charm","Spettro neutrini e interagenti CHARM CCDIS", 400, 0, 400);
 TH1D *hspectrum_nuebar_charm = new TH1D("hspectrum_nuebar_charm","Spettro antineutrini e interagenti CHARM CCDIS", 400, 0, 400);
 
 TH1D *hspectrum_numu_charm = new TH1D("hspectrum_numu_charm","Spettro neutrini mu interagenti CHARM CCDIS", 400, 0, 400);
 TH1D *hspectrum_numubar_charm = new TH1D("hspectrum_numubar_charm","Spettro antineutrini mu interagenti CHARM CCDIS", 400, 0, 400);

 TH1D *hspectrum_nue_distotal = new TH1D("hspectrum_nue_dis","Spettro neutrini e interagenti CHARM CCDIS", 400, 0, 400);
 TH1D *hspectrum_nuebar_distotal = new TH1D("hspectrum_nuebar_dis","Spettro antineutrini e interagenti CHARM CCDIS", 400, 0, 400);
 
 TH1D *hspectrum_numu_distotal = new TH1D("hspectrum_numu_dis","Spettro neutrini mu interagenti CCDIS", 400, 0, 400);
 TH1D *hspectrum_numubar_distotal = new TH1D("hspectrum_numubar_dis","Spettro antineutrini mu interagenti CCDIS", 400, 0, 400);
 
 Float_t Ninteracting = 0.;
 Double_t mass = 8183*1e+3; //mass in grams
 Double_t surface = 90.3 * 74.9; //surface in square centimetres
 //Double_t surface = 1.4e+4;

 Double_t avogadro = 6.022e+23;
 Double_t NT = mass/208. * avogadro;  
 Double_t x, y ,ysec;
 for (int n = 0; n < hfluxnutau->GetNbinsX(); n++){

  //TAU NEUTRINOS
  y = hfluxnutau->GetBinContent(n+1);
  x = hfluxnutau->GetBinLowEdge(n+1) + hfluxnutau->GetBinWidth(n+1)/2.;

  //DIS
  ysec = (hxsec_nutaup->Eval(x) + hxsec_nutaun->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_nutau->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnutau->GetBinWidth(n+1);
  hspectrum_nutau->SetBinContent(n+1,Ninteracting);
 
  //QE
  ysec = hxsec_nutau_qel_n->Eval(x); 
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnutau->GetBinWidth(n+1);
  hspectrum_nutau_qel->SetBinContent(n+1,Ninteracting);

  //MUNEUTRINOS
  y = hfluxnumu->GetBinContent(n+1);
  x = hfluxnumu->GetBinLowEdge(n+1) + hfluxnumu->GetBinWidth(n+1)/2.;
  
  //CHARMCCDIS
  ysec = (hxsec_numu_charmp->Eval(x) + hxsec_numu_charmn->Eval(x)); //summing proton and neutron cross section
  hxsec_charmtotal_numu->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnumu->GetBinWidth(n+1);
  hspectrum_numu_charm->SetBinContent(n+1,Ninteracting);

  //CCDIS
  ysec = (hxsec_numu_distotalp->Eval(x) + hxsec_numu_distotaln->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_numu->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnumu->GetBinWidth(n+1);
  hspectrum_numu_distotal->SetBinContent(n+1,Ninteracting);
 

  //TAUANTINEUTRINOS
  y = hfluxnutaubar->GetBinContent(n+1);
  x = hfluxnutaubar->GetBinLowEdge(n+1) + hfluxnutaubar->GetBinWidth(n+1)/2.;

  //DIS
  ysec = (hxsec_nutaubarp->Eval(x) + hxsec_nutaubarn->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_nutaubar->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnutaubar->GetBinWidth(n+1);
  hspectrum_nutaubar->SetBinContent(n+1,Ninteracting);
  //QE
  ysec = hxsec_nutaubar_qel_p->Eval(x); 
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnutaubar->GetBinWidth(n+1);
  hspectrum_nutaubar_qel->SetBinContent(n+1,Ninteracting);

  //MUANTINEUTRINOS
  y = hfluxnumubar->GetBinContent(n+1);
  x = hfluxnumubar->GetBinLowEdge(n+1) + hfluxnumubar->GetBinWidth(n+1)/2.;

  //CCDIS
  ysec = (hxsec_numubar_distotalp->Eval(x) + hxsec_numubar_distotaln->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_numubar->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnumubar->GetBinWidth(n+1);
  hspectrum_numubar_distotal->SetBinContent(n+1,Ninteracting);
  //CHARM CCDIS
  ysec = (hxsec_numubar_charmp->Eval(x) + hxsec_numubar_charmn->Eval(x)); //summing proton and neutron cross section
  hxsec_charmtotal_numubar->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnumubar->GetBinWidth(n+1);
  hspectrum_numubar_charm->SetBinContent(n+1,Ninteracting);

  //ELECTRON NEUTRINOS

  //ENEUTRINOS
  y = hfluxnue->GetBinContent(n+1);
  x = hfluxnue->GetBinLowEdge(n+1) + hfluxnue->GetBinWidth(n+1)/2.;

  //CCDIS
  ysec = (hxsec_nue_distotalp->Eval(x) + hxsec_nue_distotaln->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_nue->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnue->GetBinWidth(n+1);
  hspectrum_nue_distotal->SetBinContent(n+1,Ninteracting);
  
  //CHARMCCDIS
  ysec = (hxsec_nue_charmp->Eval(x) + hxsec_nue_charmn->Eval(x)); //summing proton and neutron cross section
  hxsec_charmtotal_nue->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnue->GetBinWidth(n+1);
  hspectrum_nue_charm->SetBinContent(n+1,Ninteracting);
  
  //EANTINEUTRINOS
  y = hfluxnuebar->GetBinContent(n+1);
  x = hfluxnuebar->GetBinLowEdge(n+1) + hfluxnuebar->GetBinWidth(n+1)/2.;

  //CCDIS
  ysec = (hxsec_nuebar_distotalp->Eval(x) + hxsec_nuebar_distotaln->Eval(x)); //summing proton and neutron cross section
  hxsec_distotal_nuebar->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnuebar->GetBinWidth(n+1);
  hspectrum_nuebar_distotal->SetBinContent(n+1,Ninteracting);
  
  //CHARM CCDIS
  ysec = (hxsec_nuebar_charmp->Eval(x) + hxsec_nuebar_charmn->Eval(x)); //summing proton and neutron cross section
  hxsec_charmtotal_nuebar->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnuebar->GetBinWidth(n+1);
  hspectrum_nuebar_charm->SetBinContent(n+1,Ninteracting);

 } 

 TCanvas *c2 = new TCanvas();
 hxsec_distotal_nutau->Draw();
 hxsec_distotal_nutau->SetTitle("DIS tau neutrino cross section per nucleus (#nu -> A)");
 hxsec_distotal_nutau->GetXaxis()->SetTitle("GeV");
 hxsec_distotal_nutau->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_distotal_nutaubar->SetLineColor(kRed);
 hxsec->SetLineColor(kYellow);
 hxsec_distotal_nutaubar->Draw("SAME");
 hxsec_distotal_nutau->SetFillColor(kWhite);
 hxsec_distotal_nutaubar->SetFillColor(kWhite);
 hxsec->Draw("SAME");

 TLegend *dislegend = new TLegend(0.4,0.6,0.89,0.89);
 dislegend->AddEntry("operaxsec","DIS nutau cross section formula");
 dislegend->AddEntry("dis_cc_tot_nutau", "DIS nutau cross section from GENIE");
 dislegend->AddEntry("dis_cc_tot_nutaubar", "DIS antinutau cross section from GENIE");

 dislegend->Draw("SAME");

 TCanvas *c2_1 = new TCanvas();
 hxsec_charmtotal_numu->Draw();
 hxsec_charmtotal_numu->SetTitle("CHARM CCDIS muon neutrino cross section per nucleus (#nu -> A)");
 hxsec_charmtotal_numu->GetXaxis()->SetTitle("GeV");
 hxsec_charmtotal_numu->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_charmtotal_numubar->SetLineColor(kRed);
 hxsec_charmtotal_numubar->Draw("SAME");
 hxsec_charmtotal_numu->SetFillColor(kWhite);
 hxsec_charmtotal_numubar->SetFillColor(kWhite);

 TLegend *charmmulegend = new TLegend(0.4,0.6,0.89,0.89);
 charmmulegend->AddEntry("charm_tot_numu", "CHARM CCDIS nu cross section from GENIE");
 charmmulegend->AddEntry("charm_tot_numubar", "CHARM CCDIS antinumu cross section from GENIE");
 charmmulegend->Draw("SAME");

 TCanvas *c12 = new TCanvas();
 hxsec_distotal_numu->Draw();
 hxsec_distotal_numu->SetTitle("CCDIS muon neutrino cross section per nucleus (#nu -> A)");
 hxsec_distotal_numu->GetXaxis()->SetTitle("GeV");
 hxsec_distotal_numu->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_distotal_numubar->SetLineColor(kRed);
 hxsec_distotal_numubar->Draw("SAME");
 hxsec_distotal_numu->SetFillColor(kWhite);
 hxsec_distotal_numubar->SetFillColor(kWhite);

 TLegend *distotalmulegend = new TLegend(0.4,0.6,0.89,0.89);
 distotalmulegend->AddEntry("distotal_tot_numu", "CCDIS numu cross section from GENIE");
 distotalmulegend->AddEntry("distotal_tot_numubar", "CCDIS antinumu cross section from GENIE");
 distotalmulegend->Draw("SAME");

 TCanvas *c3 = new TCanvas();
 hspectrum_nutau->Draw();
 TCanvas *c4 = new TCanvas();
 hspectrum_nutaubar->Draw();

 TCanvas *c5 = new TCanvas();
 hxsec_nutau_qel_n->Draw();
 hxsec_nutau_qel_n->SetFillColor(kWhite);
 hxsec_nutau_qel_n->SetTitle("QEL cross section tau neutrino over Pb208");
 hxsec_nutau_qel_n->GetXaxis()->SetTitle("GeV");
 hxsec_nutau_qel_n->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_nutaubar_qel_p->SetLineColor(kRed);
 hxsec_nutaubar_qel_p->SetFillColor(kWhite);
 hxsec_nutaubar_qel_p->Draw("SAME");

 TLegend *qellegend = new TLegend(0.4,0.6,0.89,0.89);
 qellegend->AddEntry("qel_cc_n", "QE numu cross section from GENIE");
 qellegend->AddEntry("qel_cc_p", "QE antinumu cross section from GENIE");
 qellegend->Draw("SAME");

 TCanvas *c2_2 = new TCanvas();
 hxsec_charmtotal_nue->Draw();
 hxsec_charmtotal_nue->SetTitle("CHARM CCDIS electron neutrino cross section per nucleus (#nu -> A)");
 hxsec_charmtotal_nue->GetXaxis()->SetTitle("GeV");
 hxsec_charmtotal_nue->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_charmtotal_nuebar->SetLineColor(kRed);
 hxsec_charmtotal_nuebar->Draw("SAME");
 hxsec_charmtotal_nue->SetFillColor(kWhite);
 hxsec_charmtotal_nuebar->SetFillColor(kWhite);

 TLegend *charmelegend = new TLegend(0.4,0.6,0.89,0.89);
 charmelegend->AddEntry("charm_tot_nue", "CHARM CCDIS nu cross section from GENIE");
 charmelegend->AddEntry("charm_tot_nuebar", "CHARM CCDIS antinumu cross section from GENIE");
 charmelegend->Draw("SAME");

 TCanvas *c11 = new TCanvas();
 hxsec_distotal_nue->Draw();
 hxsec_distotal_nue->SetTitle("CCDIS electron neutrino cross section per nucleus (#nu -> A)");
 hxsec_distotal_nue->GetXaxis()->SetTitle("GeV");
 hxsec_distotal_nue->GetYaxis()->SetTitle("10^{-38} cm^{2}");
 hxsec_distotal_nuebar->SetLineColor(kRed);
 hxsec_distotal_nuebar->Draw("SAME");
 hxsec_distotal_nue->SetFillColor(kWhite);
 hxsec_distotal_nuebar->SetFillColor(kWhite);

 TLegend *distotalelegend = new TLegend(0.4,0.6,0.89,0.89);
 distotalelegend->AddEntry("dis_tot_nue", "CCDIS nu cross section from GENIE");
 distotalelegend->AddEntry("dis_tot_nuebar", "CCDIS antinue cross section from GENIE");
 distotalelegend->Draw("SAME");


 TCanvas *c6 = new TCanvas();
 hspectrum_numu_charm->Draw();
 hspectrum_numubar_charm->Draw("SAMES");

 TCanvas *c7 = new TCanvas();
 hspectrum_nue_charm->Draw();
 hspectrum_nuebar_charm->Draw("SAMES");
 cout<<"Number of nutau DIS interactions per spill: " <<hspectrum_nutau->Integral()<<endl;
 cout<<"Number of anti nutau DIS interactions per spill: " <<hspectrum_nutaubar->Integral()<<endl;
 cout<<"Number of nutau QEL interactions per spill: " <<hspectrum_nutau_qel->Integral()<<endl;
 cout<<"Number of anti nutau QEL interactions per spill: " <<hspectrum_nutaubar_qel->Integral()<<endl;
 cout<<endl;
 cout<<"Number of nue CHARMCCDIS interactions per spill: " <<hspectrum_nue_charm->Integral()<<"over CCDIS: "<<hspectrum_nue_charm->Integral()/hspectrum_nue_distotal->Integral()<<endl;
 cout<<"Number of anti nue CHARMCCDIS interactions per spill: " <<hspectrum_nuebar_charm->Integral()<<"over CCDIS: "<<hspectrum_nuebar_charm->Integral()/hspectrum_nuebar_distotal->Integral()<<endl;;
 cout<<endl;
 cout<<"Number of nue CCDIS interactions per spill: " <<hspectrum_nue_distotal->Integral()<<endl;
 cout<<"Number of anti nue CCDIS interactions per spill: " <<hspectrum_nuebar_distotal->Integral()<<endl;
 cout<<endl;
 cout<<"Number of numu CHARMCCDIS interactions per spill: " <<hspectrum_numu_charm->Integral()<<"over CCDIS: "<<hspectrum_numu_charm->Integral()/hspectrum_numu_distotal->Integral()<<endl;
 cout<<"Number of anti numu CHARMCCDIS interactions per spill: " <<hspectrum_numubar_charm->Integral()<<"over CCDIS: "<<hspectrum_numubar_charm->Integral()/hspectrum_numubar_distotal->Integral()<<endl;
 cout<<endl;
 cout<<"Number of numu CCDIS interactions per spill: " <<hspectrum_numu_distotal->Integral()<<endl;
 cout<<"Number of anti numu CCDIS interactions per spill: " <<hspectrum_numubar_distotal->Integral()<<endl;
}





