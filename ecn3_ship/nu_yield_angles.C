//adding studies at different angles

#include "TLegend.h"
//#include "GenieGenerator.h"
#include "TGeoBBox.h"
#include <map>
//launch with root -l
//>>.L nu_yield.C
//generate_neutrinos()
const Float_t txoffset = -0.1;
void generate_neutrinos(){ //generate neutrino produced spectra according to Thomas histograms
 cout<<"chosen tx offset, propagating neutrinos "<<txoffset<<endl;

 TFile * fInputFile = TFile::Open("pythia8_Geant4_1.0_withCharm_nu.root");
 map<Int_t, TH1D*> hnu_p;
 //getting 1D spectra
 hnu_p[12] =  (TH1D*) fInputFile->Get("1012");
 hnu_p[14] =  (TH1D*) fInputFile->Get("1014");
 hnu_p[16] =  (TH1D*) fInputFile->Get("1016");
 hnu_p[-12] =  (TH1D*) fInputFile->Get("2012");
 hnu_p[-14] =  (TH1D*) fInputFile->Get("2014");
 hnu_p[-16] =  (TH1D*) fInputFile->Get("2016");

 //Float_t deltaz = 3969.; //distance between center of proton target and center of neutrino target
 Float_t deltaz = 3500.; //for now, let us assume exactly 40 m distance target-nutarget (we do not have the geometry yet)
 //Float_t deltaz = 4039.; //distance between start of proton target and center of neutrino target
 Double_t targetdx = 20.; //for geometrical acceptance requirement
 Double_t targetdy = targetdx;

 Float_t xcenter = txoffset * deltaz; //position of detector center in the x axis

 Float_t pzv;
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
 map<Int_t, TH1D*> hradiusspectrum;
 map<Int_t, Double_t> nall = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos produced, mapped per pdg
 map<Int_t, Double_t> ndet = {{12, 0.},{-12,0.},{14,0.},{-14,0.},{16,0.},{-16,0.}}; //neutrinos arrived at det

 TFile *outfile = new TFile(Form("neutrinos_detector_txoffset_%i.root",int(txoffset * 1000)),"RECREATE"); 

 hspectrumdet[12] = new TH1D("hnu_e","Spectrum electron neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-12] = new TH1D("hnu_e_bar","Spectrum electron antineutrinos arrived at detector;P[GeV/c]",400,0,400);

 hspectrumdet[14] = new TH1D("hnu_mu","Spectrum muon neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-14] = new TH1D("hnu_mu_bar","Spectrum muon antineutrinos arrived at detector;P[GeV/c]",400,0,400);

 hspectrumdet[16] = new TH1D("hnu_tau","Spectrum tau neutrinos arrived at detector;P[GeV/c]",400,0,400);
 hspectrumdet[-16] = new TH1D("hnu_tau_bar","Spectrum tau antineutrinos arrived at detector;P[GeV/c]",400,0,400);

 hradiusspectrum[12] = new TH1D("hRnu_e","Spectrum electron neutrinos arrived at detector;R[m]",200,0,20.);
 hradiusspectrum[-12] = new TH1D("hRnu_e_bar","Spectrum electron antineutrinos arrived at detector;R[m]",200,0,20.);

 hradiusspectrum[14] = new TH1D("hRnu_mu","Spectrum muon neutrinos arrived at detector;R[m]",200,0,20.);
 hradiusspectrum[-14] = new TH1D("hRnu_mu_bar","Spectrum muon antineutrinos arrived at detector;R[m]",200,0,20.);

 hradiusspectrum[16] = new TH1D("hRnu_tau","Spectrum tau neutrinos arrived at detector;R[m]",200,0,20.);
 hradiusspectrum[-16] = new TH1D("hRnu_tau_bar","Spectrum tau antineutrinos arrived at detector;R[m]",200,0,20.);

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

  if(TMath::Abs(start[0]-xcenter) < targetdx && TMath::Abs(start[1]) < targetdy){ //checking how many neutrinos are inside the detector
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
///////////////STARTING METHOD FOR ESTIMATING NEUTRINO YIELDS/////////////////////////////////////////////////////
Double_t nu_yield_general(const char* nu = "nu_mu", const char* intmode = "dis_cc", const char* charmmode = "");
void nu_yield(){ //passing neutrino types and interaction mode to nu_yield_general for estimation
  TCanvas *cspectra = new TCanvas();
  cspectra->Divide(3,3);
  cout<<"chosen tx offset, estimating yields for neutrinos "<<txoffset<<endl;
    
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
  TFile *flux = TFile::Open(Form("neutrinos_detector_txoffset_%i.root",int(txoffset * 1000)));

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
  Double_t mass = 1000*1e+3; //mass in grams
  Double_t surface = 40. * 40.; //surface in square centimetres (squared configuration)
  //Double_t surface = 1.4e+4;
   
  Double_t avogadro = 6.022e+23;
  Double_t NT = mass/A * avogadro;  
  Double_t x, y ,ysec;

  TFile *outputfile = new TFile(Form("plots_%i/results_%s_%s%s.root",int(txoffset * 1000),nu,intmode,charmmode),"RECREATE");
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

using namespace ROOT;

RVec<double> GetNuFluxes(string filename){
  //getting numbers of neutrinos in detector
  TFile * fluxfile = TFile::Open(filename.data());
  if (!fluxfile){
    cout<<"missing fluxfile check name "<<filename.data()<<endl;
  }
  const int Nnu = 6; //3 flavours, nu and antinu
  string histonames[Nnu] = {"hnu_e","hnu_mu","hnu_tau","hnu_e_bar","hnu_mu_bar","hnu_tau_bar"};
  RVec<double> nufluxes;
  for (int inu = 0; inu < Nnu; inu++){
    TH1D *hnu = (TH1D*) fluxfile->Get(histonames[inu].data());
    if (hnu) nufluxes.push_back(hnu->Integral());
    else cout<<"Missing histogram object, check name "<<histonames[inu].data()<<endl;
  }

  fluxfile->Close();
  return nufluxes;

 
}

RVec<double> GetNuYields(string foldername){
  //getting numbers of neutrinos in detector
  const int Nnu = 6; //3 flavours, nu and antinu
  string nunames[Nnu] = {"nu_e_dis_cc","nu_mu_dis_cc","nu_tau_dis_cc","nu_e_bar_dis_cc","nu_mu_bar_dis_cc","nu_tau_bar_dis_cc"};
  string histonames[Nnu] = {
    "hspectrum_nu_e_intdis_cc","hspectrum_nu_mu_intdis_cc","hspectrum_nu_tau_intdis_cc",
    "hspectrum_nu_e_bar_intdis_cc","hspectrum_nu_mu_bar_intdis_cc","hspectrum_nu_tau_bar_intdis_cc"
  };
  RVec<double> nuyields;
  for (int inu = 0; inu < Nnu; inu++){
    string prefix = "results_";
    string fileformat = ".root";
    TFile * fluxfile = TFile::Open((foldername+prefix+nunames[inu]+fileformat).data());
    if (!fluxfile){
    cout<<"missing fluxfile check name "<<(foldername+prefix+nunames[inu]+fileformat).data()<<endl;
    }

    string histoprefix = "hspectrum_";
    TH1D *hnu = (TH1D*) fluxfile->Get(histonames[inu].data());
    if (hnu) nuyields.push_back(hnu->Integral());
    else cout<<"Missing histogram object, check name "<<histonames[inu].data()<<endl;

    fluxfile->Close();
  }
  return nuyields;
}

//plot as SHiP CDS radius neutrinos
void drawAngularPlot(){
 double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
 double normship = 2e+20; //reference for five years of SND DataTaking
 const int nbins = 21;
 //listing results for different radii
 //list of folder names
 RVec<string> foldernames = {
      "plots_-100/","plots_-90/","plots_-80/","plots_-70/","plots_-60/","plots_-50/","plots_-40/","plots_-30/","plots_-20/","plots_-10/",
      "plots_0/",
      "plots_10/","plots_20/","plots_30/","plots_40/","plots_50/","plots_60/","plots_70/","plots_80/","plots_90/","plots_100/"
      };
 RVec<double> TXoffset = {-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
 //target arrays
 RVec<double> nue_ccdis,numu_ccdis,nutau_ccdis, nue_bar_ccdis, numu_bar_ccdis, nutau_bar_ccdis;
 //for each angle, I open the file and get the plots
 for (int ibin = 0; ibin < nbins; ibin++){
  RVec<double> nuyields = GetNuYields(foldernames[ibin]);
  //storing information for this angle
  nue_ccdis.push_back(nuyields[0]);
  numu_ccdis.push_back(nuyields[1]);
  nutau_ccdis.push_back(nuyields[2]);
  nue_bar_ccdis.push_back(nuyields[3]);
  numu_bar_ccdis.push_back(nuyields[4]);
  nutau_bar_ccdis.push_back(nuyields[5]);
 }
 //doing the graph. NOTE: I ADD TOGETHER NEUTRINOS AND ANTINEUTRINOS
 TGraph *gnue_ccdis = new TGraph(21, TXoffset.data(),(nue_ccdis + nue_bar_ccdis).data());
 TGraph *gnumu_ccdis = new TGraph(21, TXoffset.data(),(numu_ccdis + numu_bar_ccdis).data());
 TGraph *gnutau_ccdis = new TGraph(21, TXoffset.data(),(nutau_ccdis + nutau_bar_ccdis).data());

 gnue_ccdis->SetTitle("Electron neutrino and antineutrino;TXoffset");
 gnumu_ccdis->SetTitle("Muon neutrino and antineutrino;TXoffset");
 gnutau_ccdis->SetTitle("Tau neutrino and antineutrino;TXoffset");

 //normalizing to ship data taking
 gnue_ccdis->Scale(normship/normsim);
 gnumu_ccdis->Scale(normship/normsim);
 gnutau_ccdis->Scale(normship/normsim);

 TCanvas *cgraph = new TCanvas();
 gnumu_ccdis->SetMarkerColor(kRed);
 gnumu_ccdis->SetMarkerStyle(20);
 gnumu_ccdis->Draw("AP");
 gnumu_ccdis->GetYaxis()->SetRangeUser(1,5e+6);
 gnue_ccdis->SetMarkerColor(kBlue);
 gnue_ccdis->SetMarkerStyle(32);
 gnue_ccdis->Draw("P && SAME");
 gnutau_ccdis->SetMarkerColor(kBlack);
 gnutau_ccdis->Draw("P && SAME"); 
 gnutau_ccdis->SetMarkerStyle(21);
 cgraph->BuildLegend();
 cgraph->SetLogy();
 cgraph->Draw("g");

 gnumu_ccdis->SetTitle("CCDIS Yields with mass of 1 ton square detector at different positions");
}