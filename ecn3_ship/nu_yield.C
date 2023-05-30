#include "TLegend.h"
//#include "GenieGenerator.h"
#include "TGeoBBox.h"
#include <map>
//launch with root -l
//>>.L nu_yield.C
//generate_neutrinos()
//now with an options for logbinnins!
void generate_neutrinos(int neutrinosource = 2, bool uselogbins=false){ //generate neutrino produced spectra according to Thomas histograms

 //400,0,400
 int nbins = 400;
 double Emin = 0.;
 double Emax = 400.;
 if (uselogbins){
  nbins = 30;
  Emin = 0.;
  Emax = 3.;
 }

 cout<<"Propagating neutrino spectra with option (0 no charm, 1 only charm, 2 with charm): "<<neutrinosource<<endl;
 TString inputfilenames[3] = {"pythia8_Geant4_1.0_c_nu.root","pythia8_Geant4_charm_nu_1.0.root","pythia8_Geant4_1.0_withCharm_nu.root"};

 TFile * fInputFile = TFile::Open(inputfilenames[neutrinosource].Data());
 map<Int_t, TH1D*> hnu_p;
 //getting 1D spectra
 if(neutrinosource > 1){
  hnu_p[12] =  (TH1D*) fInputFile->Get("1012");
  hnu_p[14] =  (TH1D*) fInputFile->Get("1014");
  hnu_p[16] =  (TH1D*) fInputFile->Get("1016");
  hnu_p[-12] =  (TH1D*) fInputFile->Get("2012");
  hnu_p[-14] =  (TH1D*) fInputFile->Get("2014");
  hnu_p[-16] =  (TH1D*) fInputFile->Get("2016");
 }
 else{
  hnu_p[12] =  (TH1D*) fInputFile->Get("12");
  hnu_p[14] =  (TH1D*) fInputFile->Get("14");
  hnu_p[16] =  (TH1D*) fInputFile->Get("16");
  hnu_p[-12] =  (TH1D*) fInputFile->Get("-12");
  hnu_p[-14] =  (TH1D*) fInputFile->Get("-14");
  hnu_p[-16] =  (TH1D*) fInputFile->Get("-16");
 }

 //Float_t deltaz = 3969.; //distance between center of proton target and center of neutrino target
 Float_t deltaz = 2500.; //for now, let us assume exactly 40 m distance target-nutarget (we do not have the geometry yet)
 //Float_t deltaz = 4039.; //distance between start of proton target and center of neutrino target
 Double_t targetdx = 20.; //for geometrical acceptance requirement
 Double_t targetdy = 20.;

 Float_t pzv;
 Double_t pt = -10000.;
 Double_t start[3];

 Int_t idbase=1200;
 char ts[20];
 TH1D* pxhist[3000];//!
 TH1D* pyslice[3000][100];//!
 printf("Reading (log10(p),log10(pt)) Hists from file: %s\n",fInputFile->GetName());

//gRandom->SetSeed(0); //set 0 to make it change everytime

 for (Int_t idnu=12;idnu<17;idnu+=2){
    for (Int_t idadd=-1;idadd<2;idadd+=2){
  	  Int_t idhnu=idbase+idnu;
      if (idadd<0) idhnu+=1000; //note, we use it in all cases, even if not present in histo names. It is used to map arrays later
      if (neutrinosource > 1) sprintf(ts,"%d",idhnu); //in the with charm files, 2D histograms have names 1212, 2212
      else sprintf(ts,"2d25_%d",idnu*idadd); //in the no charm and only charm files, 2D histograms have 2d25_12, 2d25_-12 names
	    //pickup corresponding (log10(p),log10(pt)) histogram
      if (fInputFile->FindObjectAny(ts)){
           TH2F* h2tmp = (TH2F*) fInputFile->Get(ts);
           printf("HISTID=%s, Title:%s\n",ts,h2tmp->GetTitle());
	         sprintf(ts,"px_%d",idhnu);
           //make its x-projection, to later be able to convert log10(p) to its bin-number
           pxhist[idhnu]=h2tmp->ProjectionX(ts,1,-1);
           Int_t nbinx=h2tmp->GetNbinsX();
           //printf("idhnu=%d  ts=%s  nbinx=%d\n",idhnu,ts,nbinx);
	         //project all slices on the y-axis
           for (Int_t k=1;k<nbinx+1;k+=1){
	           sprintf(ts,"h%d%d",idhnu,k);
             //printf("idnu %d idhnu %d bin%d  ts=%s\n",idnu,idhnu,k,ts);
             pyslice[idhnu][k]=h2tmp->ProjectionY(ts,k,k);
	  		   }
      } //end if find object
	  } //end for loop over neutrino/antineutrinos
   } //end for loop over neutrino flavours

 //loop over histogram entries


 int neutrinopdgs[6] = {14,12,16,-14,-12,-16};
 //histogram and counters to be filled
 map<Int_t, TH1D*> hthetanu;
 map<Int_t, TH1D*> htxnu;
 map<Int_t, TH1D*> hspectrumdet;
 map<Int_t, Double_t> nall = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos produced, mapped per pdg
 map<Int_t, Double_t> ndet = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos arrived at det

 TString outputfilenames[3] = {"neutrinos_detector_nocharm.root","neutrinos_detector_charm.root","neutrinos_detector.root"};
 if (uselogbins) outputfilenames[neutrinosource] = "logbinning_"+outputfilenames[neutrinosource];
 TFile *outfile = new TFile(outputfilenames[neutrinosource],"RECREATE"); 

 hspectrumdet[12] = new TH1D("hnu_e","Spectrum electron neutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);
 hspectrumdet[-12] = new TH1D("hnu_e_bar","Spectrum electron antineutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);
 htxnu[12] = new TH1D("htx_nu_e","Electron Neutrino angular spectrum; TX",20,-0.1,0.1);
 htxnu[-12] = new TH1D("htx_nu_e_bar","Electron Antineutrino angular spectrum; TX",20,-0.1,0.1);

 hspectrumdet[14] = new TH1D("hnu_mu","Spectrum muon neutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);
 hspectrumdet[-14] = new TH1D("hnu_mu_bar","Spectrum muon antineutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);

 htxnu[14] = new TH1D("htx_nu_mu","Muon Neutrino angular spectrum; TX",20,-0.1,0.1);
 htxnu[-14] = new TH1D("htx_nu_mu_bar","Muon Antineutrino angular spectrum; TX",20,-0.1,0.1);

 hspectrumdet[16] = new TH1D("hnu_tau","Spectrum tau neutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);
 hspectrumdet[-16] = new TH1D("hnu_tau_bar","Spectrum tau antineutrinos arrived at detector;P[GeV/c]",nbins,Emin,Emax);

 htxnu[16] = new TH1D("htx_nu_tau","Tau Neutrino angular spectrum; TX",20,-0.1,0.1);
 htxnu[-16] = new TH1D("htx_nu_tau_bar","Tau Antineutrino angular spectrum; TX",20,-0.1,0.1); 
 
 //theta angles neutrinos
 hthetanu[12] = new TH1D("htheta_nu_e","Electron Neutrino angular spectrum; #theta[rad]",10,0.,0.1);
 hthetanu[-12] = new TH1D("htheta_nu_e_bar","Electron Antineutrino angular spectrum; #theta[rad]",10,0.,0.1);

 hthetanu[14] = new TH1D("htheta_nu_mu","Muon Neutrino angular spectrum; #theta[rad]",10,0.,0.1);
 hthetanu[-14] = new TH1D("htheta_nu_mu_bar","Muon Antineutrino angular spectrum; #theta[rad]",10,0.,0.1);

 hthetanu[16] = new TH1D("htheta_nu_tau","Tau Neutrino angular spectrum; #theta[rad]",10,0.,0.1);
 hthetanu[-16] = new TH1D("htheta_nu_tau_bar","Tau Antineutrino angular spectrum; #theta[rad]",10,0.,0.1); 
 

for (auto &neu:neutrinopdgs){ //start loop over neutrino flavours
 int Nentries = hnu_p[neu]->GetEntries();
 Int_t idhnu=TMath::Abs(neu)+idbase;
 if (neu<0) idhnu+=1000;
 cout<<neu<<" "<<idhnu<<endl;
 cout<<"Nentries totali: "<<Nentries<<endl;
 Double_t w = hnu_p[neu]->Integral()/hnu_p[neu]->GetEntries(); //i assume each neutrino to have the same weight (I cannot infer the original weights)
 
 for (int i = 0; i < Nentries; i++){
  pzv = hnu_p[neu]->GetRandom(); //getting p randomly from 1D histogram (pzv==pz in GENIE reference == p)
  //if (i%100000==0) cout<<i<<endl;

  // Incoming neutrino, get a random px,py
  Double_t pout[3];
  pout[2]=-1.;
  Double_t txnu=0;
  Double_t tynu=0;
  while (pout[2]<0.) {
      //Get pt of this neutrino from 2D hists.
      //Int_t idhnu=TMath::Abs(neu)+idbase;
     // if (neu<0) idhnu+=1000;
      //cout<<neu<<" "<<idhnu<<endl;
      Int_t nbinmx=pxhist[idhnu]->GetNbinsX();
      Double_t pl10=log10(pzv);
      Int_t nbx=pxhist[idhnu]->FindBin(pl10);
      //printf("idhnu %d, p %f log10(p) %f bin,binmx %d %d \n",idhnu,pzv,pl10,nbx,nbinmx);
      if (nbx<1) nbx=1;
      if (nbx>nbinmx) nbx=nbinmx;
      Double_t ptlog10=pyslice[idhnu][nbx]->GetRandom();
      //hist was filled with: log10(pt+0.01)
      pt=pow(10.,ptlog10)-0.01;
      //rotate pt in phi:
      Double_t phi=gRandom->Uniform(0.,2*TMath::Pi());
      pout[0] = cos(phi)*pt;
      pout[1] = sin(phi)*pt;
      pout[2] = pzv*pzv-pt*pt;

      if (pout[2]>=0.) {
        pout[2]=TMath::Sqrt(pout[2]);
        if (gRandom->Uniform(-1.,1.)<0.) pout[0]=-pout[0];
        if (gRandom->Uniform(-1.,1.)<0.) pout[1]=-pout[1];

        txnu=pout[0]/pout[2];
        tynu=pout[1]/pout[2];
        //cout << "Info GenieGenerator: neutrino pxyz " << pout[0] << ", " << pout[1] << ", " << pout[2] << endl;
        //printf("param %e %e %e \n",bparam,mparam[6],mparam[7]);
       }
  }
  Double_t thetanu = TMath::ASin(pt/pzv); //adding solid angle theta to the study (pt is psintheta, and pzv==pt)
  //if (i%100000==0) cout<<pout[0]<<" "<<pout[1]<<" "<<pout[2]<<" "<<pzv<<endl;
   //end px,py,pz generation, now I can follow my previous procedure for neutrino fluxes
  htxnu[neu]->Fill(txnu, w); //I want to study the spectrum of neutrinos, even if they do not enter my detector
  hthetanu[neu]->Fill(thetanu,w);

  start[0] = txnu * deltaz; //projecting produced neutrinos to neutrino detector, aggiungere x e y non cambia significativamente il risultato
  start[1] = tynu * deltaz;
  start[2] = -deltaz;

  nall[neu] +=w;

  if(TMath::Abs(start[0]) < targetdx && TMath::Abs(start[1]) < targetdy){ //checking how many neutrinos are inside the detector
    ndet[neu] += w;
    if (uselogbins) hspectrumdet[neu]->Fill(TMath::Log10(pzv),w);
    else hspectrumdet[neu]->Fill(pzv,w);
   }
  } //end main loop over neutrinos of the same flavour
 cout<<"neutrinos of flavour "<<neu<<endl;
 cout<<nall[neu]<<endl;
 cout<<ndet[neu]<<endl; //dovrebbe escludere underflow and overflow
 cout<<"Ratio arriving at det: "<<ndet[neu]/nall[neu]<<endl;
 new TCanvas();
 hspectrumdet[neu]->Draw();
 
 }//end loop over neutrino flavours
 outfile->Write();
 outfile->Close();
}
///////////////STARTING METHOD FOR ESTIMATING NEUTRINO YIELDS/////////////////////////////////////////////////////
Double_t nu_yield_general(int neutrinosource, bool uselogbins, const char* nu = "nu_mu", const char* intmode = "dis_cc", const char* charmmode = "");
void nu_yield(int neutrinosource = 2, bool uselogbins=false){ //passing neutrino types and interaction mode to nu_yield_general for estimation
  cout<<"Interacting neutrino with source option (0 no charm, 1 only charm, 2 with charm): "<<neutrinosource<<endl;
  TCanvas *cspectra = new TCanvas();
  cspectra->Divide(3,3);
    
  cout<<"Yields per nue"<<endl;
  cspectra->cd(1);
  Double_t nue_dis_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e","dis_cc");
  cspectra->cd(2);
  Double_t nue_bar_dis_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e_bar","dis_cc");
  cspectra->cd(3);
  Double_t nue_res_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e","res_cc");
  cspectra->cd(4);
  Double_t nue_bar_res_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e_bar","res_cc");
  cspectra->cd(5);
  Double_t nue_qel_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e","qel_cc");
  cspectra->cd(6);
  Double_t nue_bar_qel_cc = nu_yield_general(neutrinosource,uselogbins,"nu_e_bar","qel_cc");
  cspectra->cd(7);
  Double_t nue_dis_cc_charm =  nu_yield_general(neutrinosource,uselogbins,"nu_e","dis_cc","_charm");
  cspectra->cd(8);
  Double_t nue_bar_dis_cc_charm = nu_yield_general(neutrinosource,uselogbins,"nu_e_bar","dis_cc","_charm");
  cspectra->cd(9);
  Double_t nue_el =  nu_yield_general(neutrinosource,uselogbins,"nu_e","ve_ccncmix");
  Double_t nue_bar_el =  nu_yield_general(neutrinosource,uselogbins,"nu_e_bar","ve_ccncmix");
  cout<<endl;


  cout<<"NU: "<<" "<<"CCDIS"<<" "<<"CHARM"<<" "<<"RES"<<" "<<"QE"<<" "<<"NUEEL"<<endl;
  cout<<"NUE"<<" "<<nue_dis_cc<<" "<<nue_dis_cc_charm<<" "<<nue_res_cc<<" "<<nue_qel_cc<<" "<<nue_el<<endl;
  cout<<"ANTINUE"<<" "<<nue_bar_dis_cc<<" "<<nue_bar_dis_cc_charm<<" "<<nue_bar_res_cc<<" "<<nue_bar_qel_cc<<" "<<nue_bar_el<<endl;
  cout<<"Yields per numu"<<endl;
  TCanvas *cspectramu = new TCanvas();
  cspectramu->Divide(3,3);

  cspectramu->cd(1);
  Double_t numu_dis_cc_charm =nu_yield_general(neutrinosource,uselogbins,"nu_mu","dis_cc","_charm");
  cspectramu->cd(2);
  Double_t numu_bar_dis_cc_charm =nu_yield_general(neutrinosource,uselogbins,"nu_mu_bar","dis_cc","_charm");
  cspectramu->cd(3);
  Double_t numu_dis_cc = nu_yield_general(neutrinosource,uselogbins,"nu_mu","dis_cc");
  cspectramu->cd(4);  
  Double_t numu_bar_dis_cc =nu_yield_general(neutrinosource,uselogbins,"nu_mu_bar","dis_cc");
  cspectramu->cd(5);  
  Double_t numu_res_cc =nu_yield_general(neutrinosource,uselogbins,"nu_mu","res_cc");
  cspectramu->cd(6);
  Double_t numu_bar_res_cc =nu_yield_general(neutrinosource,uselogbins,"nu_mu_bar","res_cc");
  cspectramu->cd(7);
  Double_t numu_qel_cc =nu_yield_general(neutrinosource,uselogbins,"nu_mu","qel_cc");
  cspectramu->cd(8);
  Double_t numu_bar_qel_cc =nu_yield_general(neutrinosource,uselogbins,"nu_mu_bar","qel_cc");
  cspectramu->cd(9);
  Double_t numu_el =nu_yield_general(neutrinosource,uselogbins,"nu_mu","ve_nc");   
  Double_t numu_bar_el =nu_yield_general(neutrinosource,uselogbins,"nu_mu_bar","ve_nc");
  cout<<endl;

  cout<<"NUMU"<<" "<<numu_dis_cc<<" "<<numu_dis_cc_charm<<" "<<numu_res_cc<<" "<<numu_qel_cc<<" "<<numu_el<<endl;
  cout<<"ANTINUMU"<<" "<<numu_bar_dis_cc<<" "<<numu_bar_dis_cc_charm<<" "<<numu_bar_res_cc<<" "<<numu_bar_qel_cc<<" "<<numu_bar_el<<endl;
  
  TCanvas *cspectratau = new TCanvas();
  cspectratau->Divide(3,3);
  cout<<"Yields per nutau"<<endl;
  cspectratau->cd(1);
  Double_t nutau_dis_cc_charm =nu_yield_general(neutrinosource,uselogbins,"nu_tau","dis_cc","_charm");
  cspectratau->cd(2);
  Double_t nutau_bar_dis_cc_charm =nu_yield_general(neutrinosource,uselogbins,"nu_tau_bar","dis_cc","_charm");
  cspectratau->cd(3);
  Double_t nutau_dis_cc = nu_yield_general(neutrinosource,uselogbins,"nu_tau","dis_cc");
  cspectratau->cd(4);
  Double_t nutau_bar_dis_cc =nu_yield_general(neutrinosource,uselogbins,"nu_tau_bar","dis_cc");
  cspectratau->cd(5);
  Double_t nutau_qel_cc =nu_yield_general(neutrinosource,uselogbins,"nu_tau","qel_cc");
  cspectratau->cd(6);
  Double_t nutau_bar_qel_cc =nu_yield_general(neutrinosource,uselogbins,"nu_tau_bar","qel_cc");
  cspectratau->cd(7);  
  Double_t nutau_res_cc =nu_yield_general(neutrinosource,uselogbins,"nu_tau","res_cc");
  cspectratau->cd(8);
  Double_t nutau_bar_res_cc =nu_yield_general(neutrinosource,uselogbins,"nu_tau_bar","res_cc");
  cspectratau->cd(9);
  Double_t nutau_el =nu_yield_general(neutrinosource,uselogbins,"nu_tau","ve_nc");   
  Double_t nutau_bar_el =nu_yield_general(neutrinosource,uselogbins,"nu_tau_bar","ve_nc");

  cout<<"NUTAU"<<" "<<nutau_dis_cc<<" "<<nutau_dis_cc_charm<<" "<<nutau_res_cc<<" "<<nutau_qel_cc<<" "<<nutau_el<<endl;
  cout<<"ANTINUTAU"<<" "<<nutau_bar_dis_cc<<" "<<nutau_bar_dis_cc_charm<<" "<<nutau_bar_res_cc<<" "<<nutau_bar_qel_cc<<" "<<nutau_bar_el<<endl;
  
}

void BinLogX(TH1 *h)
{
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();
   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];
   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
}

//general layout with FORM to estimate number of neutrino interactions 
Double_t nu_yield_general(int neutrinosource, bool uselogbins, const char* nu = "nu_mu", const char* intmode = "dis_cc", const char* charmmode = ""){     
  TFile *xsec = TFile::Open("nu_xsec_TungstenSHIP.root");
  //TFile *xsec = TFile::Open("Nu_xsec_full.root"); //normal splines are cut at 350 GeV
  TFile *flux = NULL;
  TString fluxfilenames[3] = {"neutrinos_detector_nocharm.root","neutrinos_detector_charm.root","neutrinos_detector.root"};
  if (uselogbins) fluxfilenames[neutrinosource] = "logbinning_"+fluxfilenames[neutrinosource];
  flux = TFile::Open(fluxfilenames[neutrinosource]);

  //const char* nu = "nu_mu"; //neutrino type
  //const char* intmode = "dis_cc"; //interaction mode;

  //  const char* charmmode = "_charm"; //write _charm for CHARMCCDIS
  
  TH1D *hfluxnu = (TH1D*) flux->Get(Form("h%s",nu));

  bool nullp = false, nulln = false, nullother = false; //for DIS both are possible, but not for other interactions
  //TGraph *hxsec_p = (TGraph*) xsec->Get(Form("%s_Pb208/%s_p%s",nu,intmode,charmmode));   
  //TGraph *hxsec_n = (TGraph*) xsec->Get(Form("%s_Pb208/%s_n%s",nu,intmode,charmmode));
  //TGraph *hxsec_other = (TGraph*) xsec->Get(Form("%s_Pb208/%s%s",nu,intmode,charmmode));
  TGraph *hxsec_p = (TGraph*) xsec->Get(Form("%s_W184/%s_p%s",nu,intmode,charmmode));   
  TGraph *hxsec_n = (TGraph*) xsec->Get(Form("%s_W184/%s_n%s",nu,intmode,charmmode));
  TGraph *hxsec_other = (TGraph*) xsec->Get(Form("%s_W184/%s%s",nu,intmode,charmmode));
  if (hxsec_p == NULL) nullp = true;
  if (hxsec_n == NULL) nulln = true;
  if (hxsec_other == NULL) nullother = true;
  
  //const Int_t A = 208; //initial approximation, only pb208 considered
  // const Int_t Z = 82;

  const Int_t A = 184; //initial approximation, only pb208 considered
  const Int_t Z = 74;

  Float_t Ninteracting = 0.;
  //Double_t mass = 8183*1e+3; //mass in grams
  //Double_t mass = 8352*1e+3; //mass in grams
  Double_t mass = 3080*1e+3; //mass in grams
  Double_t surface = 40. * 40.; //surface in square centimetres (squared configuration)
  //Double_t surface = 1.4e+4;
   
  Double_t avogadro = 6.022e+23;
  Double_t NT = mass/A * avogadro;  
  Double_t x, y ,ysec;

  TString plotsfolder = Form("plots_%i/results_%s_%s%s.root",neutrinosource, nu,intmode,charmmode);
  if (uselogbins) plotsfolder = "logbinning_"+plotsfolder;

  TFile *outputfile = new TFile(plotsfolder.Data(),"RECREATE");
  TGraph *hxsec_total = new TGraph(); //summing protons and neutrons when both are present
  int nbins = hfluxnu->GetNbinsX();
  double Emin = hfluxnu->GetXaxis()->GetXmin();
  double Emax = hfluxnu->GetXaxis()->GetXmax();

  //cout<<"Reading histogram with nbins "<<nbins<<" from minimum energy "<<Emin<<" to maximum energy "<<Emax<<endl;

  TH1D *hspectrum_int = new TH1D(Form("hspectrum_%s_int%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), nbins, Emin, Emax);

  for (int n = 0; n < hfluxnu->GetNbinsX(); n++){
  
  y = hfluxnu->GetBinContent(n+1);
  //if log10 was used to fill the histogram, we have to get back to their value in base 10
  if(uselogbins) x = pow(10,hfluxnu->GetXaxis()->GetBinCenter(n+1)); 
  else x = hfluxnu->GetXaxis()->GetBinCenter(n+1);


  
  if (nullp && nulln && nullother) cout<<Form("ERROR: no cross sections are present for neutrino %s in interaction %s%s",nu,intmode,charmmode)<<endl;
  if (nullp && nulln) ysec = hxsec_other->Eval(x);
  else if (nullp) ysec = hxsec_n->Eval(x);
  else if (nulln) ysec = hxsec_p->Eval(x);
  else ysec = (hxsec_p->Eval(x) + hxsec_n->Eval(x)); //summing proton and neutron cross section

  hxsec_total->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface;  
  //hspectrum_int->GetXaxis()->SetRange(1,nbins);
  hspectrum_int->SetBinContent(n+1,Ninteracting);
  }
  Double_t counting = hfluxnu->Integral();

  //Double_t yield =  hspectrum_int->Integral(Emin,Emax);
  Double_t yield =  hspectrum_int->Integral();
  if(uselogbins)BinLogX(hspectrum_int);
  //cout<<Form("numero neutrini %s per spill interagenti in mode %s%s",nu,intmode,charmmode)<<" "<<hspectrum_int->Integral()<<endl;
  hspectrum_int->Draw();
  hspectrum_int->GetXaxis()->SetTitle("GeV/c");
  // hxsec_total->Draw();
  hxsec_total->SetTitle(Form("cross section from GENIE for %s in interaction mode %s%s",nu,intmode,charmmode));
  outputfile->Write();
  hxsec_total->Write("xsec");
  outputfile->Close();
  return yield;
}

void drawInteractingSpectra(){ //drawing the spectra of interacting CCDIS
  const double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
  const double normship = 2e+20; //reference for one year of SND DataTaking

  //Drawing stored histograms for neutrinos (/home/utente/Simulations/nuyield_shipecn3/25m/)
  TFile *nuefile = TFile::Open("plots_2/results_nu_e_dis_cc.root");
  TFile *nuebarfile = TFile::Open("plots_2/results_nu_e_bar_dis_cc.root");

  TFile *numufile = TFile::Open("plots_2/results_nu_mu_dis_cc.root");
  TFile *numubarfile = TFile::Open("plots_2/results_nu_mu_bar_dis_cc.root");

  TFile *nutaufile = TFile::Open("plots_2/results_nu_tau_dis_cc.root");
  TFile *nutaubarfile = TFile::Open("plots_2/results_nu_tau_bar_dis_cc.root");

  //getting histograms
  TH1D *hnue = (TH1D*) nuefile->Get("hspectrum_nu_e_intdis_cc");
  TH1D *hnuebar = (TH1D*) nuebarfile->Get("hspectrum_nu_e_bar_intdis_cc");

  TH1D *hnumu = (TH1D*) numufile->Get("hspectrum_nu_mu_intdis_cc");
  TH1D *hnumubar = (TH1D*) numubarfile->Get("hspectrum_nu_mu_bar_intdis_cc");

  TH1D *hnutau = (TH1D*) nutaufile->Get("hspectrum_nu_tau_intdis_cc");
  TH1D *hnutaubar = (TH1D*) nutaubarfile->Get("hspectrum_nu_tau_bar_intdis_cc");

  //naming histogram
  hnue->SetTitle("Electron neutrino");
  hnumu->SetTitle("Muon neutrino");
  hnutau->SetTitle("Tau neutrino");

  hnuebar->SetTitle("Electron anti-neutrino");
  hnumubar->SetTitle("Muon anti-neutrino");
  hnutaubar->SetTitle("Tau anti-neutrino");

  //normalizing histograms to the data taking
  
  hnue->Scale(normship/normsim);
  hnuebar->Scale(normship/normsim);

  hnumu->Scale(normship/normsim);
  hnumubar->Scale(normship/normsim);

  hnutau->Scale(normship/normsim);
  hnutaubar->Scale(normship/normsim);

  cout<<"adding together neutrinos and antineutrinos of the same flavour"<<endl;
  hnue->Add(hnuebar);
  hnumu->Add(hnumubar);
  hnutau->Add(hnutaubar);

  //drawing and building legend
  TCanvas *cnuspectra = new TCanvas();
  hnumu->SetLineColor(kRed);
  hnutau->SetLineColor(kGreen);
  
  hnuebar->SetLineColor(kBlack);
  hnumubar->SetLineColor(kMagenta);
  hnutaubar->SetLineColor(kYellow);

  hnumu->Draw("histo");
  hnue->Draw("histo&&SAMES");
  hnutau->Draw("histo&&SAMES");

  //hnuebar->Draw("histo&&SAMES");
  //hnumubar->Draw("histo&&SAMES");
  //hnutaubar->Draw("histo&&SAMES");
  cnuspectra->BuildLegend();
  hnumu->SetTitle("Interacting spectra");

}

void drawSpectra(){
  const double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
  const double normship = 2e+20; //reference for five years of SND DataTaking
  //Drawing stored histograms for neutrinos (/home/utente/Simulations/nuyield_shipecn3/35m/)
  TFile *inputfile = TFile::Open("neutrinos_detector.root");
  //Getting plots
  TH1D *htx_nu_e = (TH1D*) inputfile->Get("htx_nu_e");
  TH1D *htx_nu_e_bar = (TH1D*) inputfile->Get("htx_nu_e_bar");  

  TH1D *htx_nu_mu = (TH1D*) inputfile->Get("htx_nu_mu");
  TH1D *htx_nu_mu_bar = (TH1D*) inputfile->Get("htx_nu_mu_bar");

  TH1D *htx_nu_tau = (TH1D*) inputfile->Get("htx_nu_tau");
  TH1D *htx_nu_tau_bar = (TH1D*) inputfile->Get("htx_nu_tau_bar");

  //Theta plots
  TH1D *htheta_nu_e = (TH1D*) inputfile->Get("htheta_nu_e");
  TH1D *htheta_nu_e_bar = (TH1D*) inputfile->Get("htheta_nu_e_bar");  

  TH1D *htheta_nu_mu = (TH1D*) inputfile->Get("htheta_nu_mu");
  TH1D *htheta_nu_mu_bar = (TH1D*) inputfile->Get("htheta_nu_mu_bar");

  TH1D *htheta_nu_tau = (TH1D*) inputfile->Get("htheta_nu_tau");
  TH1D *htheta_nu_tau_bar = (TH1D*) inputfile->Get("htheta_nu_tau_bar");


  //summing neutrinos and antineutrinos together
  htx_nu_e->Add(htx_nu_e_bar);
  htx_nu_mu->Add(htx_nu_mu_bar);
  htx_nu_tau->Add(htx_nu_tau_bar);

  htheta_nu_e->Add(htheta_nu_e_bar);
  htheta_nu_mu->Add(htheta_nu_mu_bar);
  htheta_nu_tau->Add(htheta_nu_tau_bar);

  //normalizing to five iyears of data taking
  htx_nu_e->Scale(normship/normsim);
  htx_nu_mu->Scale(normship/normsim);
  htx_nu_tau->Scale(normship/normsim);

  htheta_nu_e->Scale(normship/normsim);
  htheta_nu_mu->Scale(normship/normsim);
  htheta_nu_tau->Scale(normship/normsim);

  //drawing them together
  TCanvas *c = new TCanvas();
  htx_nu_mu->SetLineColor(kRed);
  htx_nu_mu->Draw();
  htx_nu_e->SetLineColor(kBlue);
  htx_nu_e->Draw("SAMES");
  htx_nu_tau->SetLineColor(kBlack);
  htx_nu_tau->Draw("SAMES");
  c->BuildLegend();

  TCanvas *ctheta = new TCanvas();
  htheta_nu_mu->SetLineColor(kRed);
  htheta_nu_mu->Draw();
  htheta_nu_e->SetLineColor(kBlue);
  htheta_nu_e->Draw("SAMES");
  htheta_nu_tau->SetLineColor(kBlack);
  htheta_nu_tau->Draw("SAMES");
  ctheta->BuildLegend();





}
