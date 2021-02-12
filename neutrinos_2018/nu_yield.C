#include "TLegend.h"
#include "GenieGenerator.h"
#include "TGeoBBox.h"
#include <map>
//when you discover the difficulty of copying from the masters. Trying to replicate Annarita's studies on neutrino fluxes

void generate_neutrinos(){ //generate neutrino produced spectra according to Thomas histograms

 TFile * fInputFile = TFile::Open("pythia8_Geant4_1.0_withCharm_nu.root");
 map<Int_t, TH1D*> hnu_p;
 //getting 1D spectra
 hnu_p[12] =  (TH1D*) fInputFile->Get("1012");
 hnu_p[14] =  (TH1D*) fInputFile->Get("1014");
 hnu_p[16] =  (TH1D*) fInputFile->Get("1016");
 hnu_p[-12] =  (TH1D*) fInputFile->Get("2012");
 hnu_p[-14] =  (TH1D*) fInputFile->Get("2014");
 hnu_p[-16] =  (TH1D*) fInputFile->Get("2016");

 Float_t deltaz = 3969.; //distance between center of proton target and start of neutrino target
 Double_t targetdx = 40., targetdy = 40.; //for geometrical acceptance requirement

 Float_t pzv;
 Double_t start[3];

 Int_t idbase=1200;
 char ts[20];
 TH1D* pxhist[3000];//!
 TH1D* pyslice[3000][100];//!
 printf("Reading (log10(p),log10(pt)) Hists from file: %s\n",fInputFile->GetName());
 for (Int_t idnu=12;idnu<17;idnu+=2){
    for (Int_t idadd=-1;idadd<2;idadd+=2){
  	  Int_t idhnu=idbase+idnu;
      if (idadd<0) idhnu+=1000;
	    sprintf(ts,"%d",idhnu);
	    //pickup corresponding (log10(p),log10(pt)) histogram
      if (fInputFile->FindObjectAny(ts)){
           TH2F* h2tmp = (TH2F*) fInputFile->Get(ts);
           printf("HISTID=%d, Title:%s\n",idhnu,h2tmp->GetTitle());
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
 map<Int_t, TH1D*> hspectrumdet;
 map<Int_t, Double_t> nall = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos produced, mapped per pdg
 map<Int_t, Double_t> ndet = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos arrived at det

 TFile *outfile = new TFile("neutrinos_detector.root","RECREATE"); 

 hspectrumdet[12] = new TH1D("hnu_e","Spectrum electron neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-12] = new TH1D("hnu_e_bar","Spectrum electron antineutrinos arrived at detector;P[GeV/c]",400,0,400);

 hspectrumdet[14] = new TH1D("hnu_mu","Spectrum muon neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-14] = new TH1D("hnu_mu_bar","Spectrum muon antineutrinos arrived at detector;P[GeV/c]",400,0,400);

 hspectrumdet[16] = new TH1D("hnu_tau","Spectrum tau neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-16] = new TH1D("hnu_tau_bar","Spectrum tau antineutrinos arrived at detector;P[GeV/c]",400,0,400);

for (auto &neu:neutrinopdgs){ //start loop over neutrino flavours
 int Nentries = hnu_p[neu]->GetEntries();
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
      Int_t idhnu=TMath::Abs(neu)+idbase;
      if (neu<0) idhnu+=1000;
      Int_t nbinmx=pxhist[idhnu]->GetNbinsX();
      Double_t pl10=log10(pzv);
      Int_t nbx=pxhist[idhnu]->FindBin(pl10);
      //printf("idhnu %d, p %f log10(p) %f bin,binmx %d %d \n",idhnu,pzv,pl10,nbx,nbinmx);
      if (nbx<1) nbx=1;
      if (nbx>nbinmx) nbx=nbinmx;
      Double_t ptlog10=pyslice[idhnu][nbx]->GetRandom();
      //hist was filled with: log10(pt+0.01)
      Double_t pt=pow(10.,ptlog10)-0.01;
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
  //if (i%100000==0) cout<<pout[0]<<" "<<pout[1]<<" "<<pout[2]<<" "<<pzv<<endl;
   //end px,py,pz generation, now I can follow my previous procedure for neutrino fluxes

  start[0] = txnu * deltaz; //projecting produced neutrinos to neutrino detector, aggiungere x e y non cambia significativamente il risultato
  start[1] = tynu * deltaz;
  start[2] = -deltaz;

  nall[neu] +=w;

  if(TMath::Abs(start[0]) < targetdx && TMath::Abs(start[1]) < targetdy){ //checking how many neutrinos are inside the detector
    ndet[neu] += w;
    hspectrumdet[neu]->Fill(pzv,w);
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

void neutrino_fluxes(){ //projecting neutrino fluxes to the target
 TFile * f = TFile::Open("/home/antonio/SHIPBuild/pythia8_Geant4-withCharm_onlyNeutrinos.root");
 TTree *tree = (TTree*) f->Get("pythia8-Geant4");
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("geofile_full.conical.Genie-TGeant4.root");
 Double_t targetZ = 2* ((TGeoBBox*) (gGeoManager->GetVolume("volTarget")->GetShape()))->GetDZ();

 //Double_t targetZ = 320.83; 

 cout<<"Spessore bersaglio neutrini da fairship: "<<targetZ<<endl; 

 Double_t bparam;
 Double_t mparam[10];

 GenieGenerator *fMaterialInvestigator = new GenieGenerator();
 Double_t start[3], end[3]; //start and end of neutrino positions at z of neutrino target
 const Int_t nneutrinos = tree->GetEntries();
 //getting the branches

 Float_t id, px, py, pz,w,x,y;

 tree->SetBranchAddress("id",&id);
 tree->SetBranchAddress("x",&x);
 tree->SetBranchAddress("y",&y);
 tree->SetBranchAddress("px",&px);
 tree->SetBranchAddress("py",&py);
 tree->SetBranchAddress("pz",&pz);
 tree->SetBranchAddress("w",&w);
 
 Float_t deltaz = 3969.; //distance between center of proton target and start of neutrino target
 //Float_t deltaz = 5500.;
 Float_t xfin, yfin, tanx,tany; //angles and positions after projections

 Float_t nnudet = 0., nallneutrinos = 0.;
 map<Int_t, Float_t> nall = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos produced, mapped per pdg
 map<Int_t, Float_t> ndet = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos arrived at det

 map<Int_t, TH1D*> hspectrumdet;

 TFile *outfile = new TFile("neutrinos_detector_upto1GeV.root","RECREATE"); //spectra arrived in detector are saved here
 //nutau
 hspectrumdet[16] = new TH1D("hnu_tau","Spectrum tau neutrinos arrived at detector",1000,0,1);
 hspectrumdet[-16] = new TH1D("hnu_tau_bar","Spectrum tau neutrinos arrived at detector",1000,0,1);
 TH2D *hxy_nutau_arrived = new TH2D("hxy", "XY distribution of tau neutrinos at detector",2000,-1000,1000,2000,-1000,1000);
 TH2D *hxy_nutau_arrived_det = new TH2D("hxy_det", "XY distribution of tau neutrinos at detector",90,-45,45,76,-38,38);

 TH1D *hnutaubar_weight = new TH1D("hnutaubar_weight", "Mean density per length transversed in target region",400,0,400);

 //numu
 hspectrumdet[14] = new TH1D("hnu_mu","Spectrum muon neutrinos arrived at detector",1000,0,1);
 hspectrumdet[-14] = new TH1D("hnu_mu_bar","Spectrum muon antineutrinos arrived at detector",1000,0,1);

 hspectrumdet[12] = new TH1D("hnu_e","Spectrum electron neutrinos arrived at detector",1000,0,1);
 hspectrumdet[-12] = new TH1D("hnu_e_bar","Spectrum electron antineutrinos arrived at detector",1000,0,1);

 Double_t targetdx = 45.15, targetdy = 37.45;
 //Double_t targetdx = 55, targetdy = 55;
 cout<<"N NEUTRINOS: "<<nneutrinos<<endl;
 for (int i = 0; i < nneutrinos; i++){

 tree->GetEntry(i);
 if (i%100000 == 0) cout<<i<<endl;
 tanx = px/pz;
 tany = py/pz;

 Double_t momentum = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

 start[0] = tanx * deltaz; //projecting produced neutrinos to neutrino detector, aggiungere x e y non cambia significativamente il risultato
 start[1] = tany * deltaz;
 start[2] = -3259;

 end[0] = tanx * (deltaz + targetZ);
 end[1] = tany * (deltaz + targetZ);
 end[2] = start[2] + targetZ;

 nallneutrinos += w;

 nall[id] +=w;

 if(TMath::Abs(start[0]) < targetdx && TMath::Abs(start[1]) < targetdy){ //checking how many neutrinos are inside the detector
  nnudet += w;
  ndet[id] += w;
  hspectrumdet[id]->Fill(momentum,w);
  bparam = fMaterialInvestigator->MeanMaterialBudget(start, end, mparam);
  if (id == -16) hnutaubar_weight->Fill(mparam[0]/mparam[4]);
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

 TCanvas *c4 = new TCanvas();
 hxy_nutau_arrived->Draw("COLZ");
 TCanvas *c5 = new TCanvas();
 hxy_nutau_arrived_det->Draw("COLZ");

 outfile->Write();
 outfile->Close();
}

///////////////STARTING METHOD FOR ESTIMATING NEUTRINO YIELDS/////////////////////////////////////////////////////
Double_t nu_yield_general(const char* nu = "nu_mu", const char* intmode = "dis_cc", const char* charmmode = "");
void nu_yield(){ //passing neutrino types and interaction mode to nu_yield_general for estimation
  TCanvas *cspectra = new TCanvas();
  cspectra->Divide(3,3);
    
  cout<<"Yields per nue"<<endl;
  cspectra->cd(1);
  Double_t nue_dis_cc = nu_yield_general("nu_e","dis_cc");
  cspectra->cd(2);
  Double_t nue_bar_dis_cc = nu_yield_general("nu_e_bar","dis_cc");
  cspectra->cd(3);
  Double_t nue_res_cc = nu_yield_general("nu_e","res_cc");
  cspectra->cd(4);
  Double_t nue_bar_res_cc = nu_yield_general("nu_e_bar","res_cc");
  cspectra->cd(5);
  Double_t nue_qel_cc = nu_yield_general("nu_e","qel_cc");
  cspectra->cd(6);
  Double_t nue_bar_qel_cc = nu_yield_general("nu_e_bar","qel_cc");
  cspectra->cd(7);
  Double_t nue_dis_cc_charm =  nu_yield_general("nu_e","dis_cc","_charm");
  cspectra->cd(8);
  Double_t nue_bar_dis_cc_charm = nu_yield_general("nu_e_bar","dis_cc","_charm");
  cspectra->cd(9);
  Double_t nue_el =  nu_yield_general("nu_e","ve_ccncmix");
  Double_t nue_bar_el =  nu_yield_general("nu_e_bar","ve_ccncmix");
  cout<<endl;


  cout<<"NU: "<<" "<<"CCDIS"<<" "<<"CHARM"<<" "<<"RES"<<" "<<"QE"<<" "<<"NUEEL"<<endl;
  cout<<"NUE"<<" "<<nue_dis_cc<<" "<<nue_dis_cc_charm<<" "<<nue_res_cc<<" "<<nue_qel_cc<<" "<<nue_el<<endl;
  cout<<"ANTINUE"<<" "<<nue_bar_dis_cc<<" "<<nue_bar_dis_cc_charm<<" "<<nue_bar_res_cc<<" "<<nue_bar_qel_cc<<" "<<nue_bar_el<<endl;
  cout<<"Yields per numu"<<endl;
  TCanvas *cspectramu = new TCanvas();
  cspectramu->Divide(3,3);

  cspectramu->cd(1);
  Double_t numu_dis_cc_charm =nu_yield_general("nu_mu","dis_cc","_charm");
  cspectramu->cd(2);
  Double_t numu_bar_dis_cc_charm =nu_yield_general("nu_mu_bar","dis_cc","_charm");
  cspectramu->cd(3);
  Double_t numu_dis_cc = nu_yield_general("nu_mu","dis_cc");
  cspectramu->cd(4);  
  Double_t numu_bar_dis_cc =nu_yield_general("nu_mu_bar","dis_cc");
  cspectramu->cd(5);  
  Double_t numu_res_cc =nu_yield_general("nu_mu","res_cc");
  cspectramu->cd(6);
  Double_t numu_bar_res_cc =nu_yield_general("nu_mu_bar","res_cc");
  cspectramu->cd(7);
  Double_t numu_qel_cc =nu_yield_general("nu_mu","qel_cc");
  cspectramu->cd(8);
  Double_t numu_bar_qel_cc =nu_yield_general("nu_mu_bar","qel_cc");
  cspectramu->cd(9);
  Double_t numu_el =nu_yield_general("nu_mu","ve_nc");   
  Double_t numu_bar_el =nu_yield_general("nu_mu_bar","ve_nc");
  cout<<endl;

  cout<<"NUMU"<<" "<<numu_dis_cc<<" "<<numu_dis_cc_charm<<" "<<numu_res_cc<<" "<<numu_qel_cc<<" "<<numu_el<<endl;
  cout<<"ANTINUMU"<<" "<<numu_bar_dis_cc<<" "<<numu_bar_dis_cc_charm<<" "<<numu_bar_res_cc<<" "<<numu_bar_qel_cc<<" "<<numu_bar_el<<endl;
  
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
  TFile *xsec = TFile::Open("nu_xsec_TungstenSHIP.root");
  //TFile *xsec = TFile::Open("Nu_xsec_full.root"); //normal splines are cut at 350 GeV
  TFile *flux = TFile::Open("neutrinos_detector.root");

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
  Double_t mass = 8352*1e+3; //mass in grams
  Double_t surface = 80. * 80.; //surface in square centimetres (squared configuration)
  //Double_t surface = 1.4e+4;
  
  Double_t avogadro = 6.022e+23;
  Double_t NT = mass/A * avogadro;  
  Double_t x, y ,ysec;

  TFile *outputfile = new TFile(Form("plots/results_%s_%s%s.root",nu,intmode,charmmode),"RECREATE");
  TGraph *hxsec_total = new TGraph(); //summing protons and neutrons when both are present
  TH1D *hspectrum_int = new TH1D(Form("hspectrum_%s_int%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), 400, 0, 400); 
  TH1D *htest = new TH1D(Form("htest_%s_int%s%s",nu,intmode,charmmode),Form("Spettro neutrini %s interagenti in %s%s",nu,intmode,charmmode), 400, 0, 400); 

  for (int n = 0; n < hfluxnu->GetNbinsX(); n++){
  
  y = hfluxnu->GetBinContent(n+1);
  x = hfluxnu->GetBinLowEdge(n+1) + hfluxnu->GetBinWidth(n+1)/2.;
  
  if (nullp && nulln && nullother) cout<<Form("ERROR: no cross sections are present for neutrino %s in interaction %s%s",nu,intmode,charmmode)<<endl;
  if (nullp && nulln) ysec = hxsec_other->Eval(x);
  else if (nullp) ysec = hxsec_n->Eval(x);
  else if (nulln) ysec = hxsec_p->Eval(x);
  else ysec = (hxsec_p->Eval(x) + hxsec_n->Eval(x)); //summing proton and neutron cross section

  hxsec_total->SetPoint(n+1,x,ysec);
  Ninteracting = ysec *1e-38 * NT * y/surface * hfluxnu->GetBinWidth(n+1);  
  htest->SetBinContent(n+1,hfluxnu->GetBinWidth(n+1)*y);
  hspectrum_int->SetBinContent(n+1,Ninteracting);

  }
  Double_t counting = hfluxnu->Integral();

  Double_t yield =  hspectrum_int->Integral(0,400);
  //cout<<Form("numero neutrini %s per spill interagenti in mode %s%s",nu,intmode,charmmode)<<" "<<hspectrum_int->Integral()<<endl;
  hspectrum_int->Draw();
  hspectrum_int->GetXaxis()->SetTitle("GeV/c");
  // hxsec_total->Draw();
  hxsec_total->SetTitle(Form("cross section from GENIE for %s in interaction mode %s%s",nu,intmode,charmmode));
  outputfile->Write();
  hxsec_total->Write("xsec");
  outputfile->Close();
  //cout<<"TEST "<<htest->GetEntries()<<endl;
  //Double_t countingtest = htest->Integral();
  //cout<<"Test conteggio "<<counting<<"integrale manuale "<<countingtest<<endl;
  return yield;
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



