//reading numu and numubar simulation, AFTER they have cut within target acceptance with cutsim.py
double GetParticleCharge (int pdgcode, TDatabasePDG *pdg){
  //from PDG, get charge
  double charge = 0.;
  if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
  else if (pdgcode > 1e+8) charge = 1.; //test storing heavy nuclei
  return charge;
}
using namespace ROOT;
//in this configuration, we have two tracking stations upstream of the magnet, two downstream. Bending Angle measurement
RVec<Double_t> BendingSpectrometer(Double_t *XSpectrometer, Double_t *YSpectrometer, Double_t * ZSpectrometer){
    //0 and 1: upstream stations;
    Double_t TXup = (XSpectrometer[1] - XSpectrometer[0])/(ZSpectrometer[1] - ZSpectrometer[0]);
    Double_t TYup = (YSpectrometer[1] - YSpectrometer[0])/(ZSpectrometer[1] - ZSpectrometer[0]);
    //2 and 3: downstream stations;
    Double_t TXdown = (XSpectrometer[3] - XSpectrometer[2])/(ZSpectrometer[3] - ZSpectrometer[2]);
    Double_t TYdown = (YSpectrometer[3] - YSpectrometer[2])/(ZSpectrometer[3] - ZSpectrometer[2]);
  
    RVec<Double_t> DeltaT = {-99999.,-99999.};

    DeltaT[0] = TXdown - TXup;
    DeltaT[1] = TYdown - TYup;

    return DeltaT;
}
//in this configuration, we have two tracking stations outside of the magnet, two inside. Sagitta measurement
RVec<Double_t> SagittaSpectrometer(Double_t *XSpectrometer, Double_t *YSpectrometer, Double_t *ZSpectrometer){
  //0 and 4: upstream stations;
  Double_t MeanXout = (XSpectrometer[0] + XSpectrometer[3])/2.;
  Double_t MeanYout = (YSpectrometer[0] + YSpectrometer[3])/2.;

  Double_t MeanXin = (XSpectrometer[1] + XSpectrometer[2])/2.;
  Double_t MeanYin = (YSpectrometer[1] + YSpectrometer[2])/2.;

  RVec<Double_t> DeltaT = {-99999.,-99999.};

  DeltaT[0] = MeanXin - MeanXout;
  DeltaT[1] = MeanYin - MeanYout;

  return DeltaT;
}
//verify muon charge separation from angular smearing along y
void muonspectrometer_testsimulation(){
 TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
 ROOT::RVec<int> charmpdglist = {421,411,431,4122,4232,4132,4332}; //to study CharmCCDIS

 const int eventesclusion = 2; //0 only odd, 1 only even, 2 all;
 const float dx_acceptance = 60.; //how many we lose by reducing our station size? (this is HALF SIZE)
 const float dy_acceptance = 60.; 

 TString prefix("root:://eosuser.cern.ch/");//for ROOTXD
 //TString simpath_mu("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/numu_CCDIS_2023_03_04_targetmovedupstream_spectromag/");
 //TString simpath_mu_bar("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/numu_bar_CCDIS_2023_03_04_targetmovedupstream_spectromag/");
 //TString simpath_mu_bar("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_04_14_numu_bar_CHARMCCDIS_spectro_1_2T/");
 TString simpath_mu("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_05_27_numu_CCDIS_spectrosagitta/");
 TRandom3 *randomgen  = new TRandom3(); //0; seed changes everytime, no seed, default to 4357 

 const double spectro_posres = 100.*1e-4; //100 micron, for smearing 
 
 TChain *simchain = new TChain("cbmsim");
 simchain->Add((prefix+simpath_mu+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 //simchain->Add((prefix+simpath_mu_bar+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 
 TTreeReader reader(simchain);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 TProfile *profmuE_nuE = new TProfile("profmuE_nuE","Profile muon energy vs neutrino energy;E_mu[GeV];E_mu/E_nu",400,0,400,0,1);

 TH1D * hothernudauP = new TH1D("hothernudauP","momentum of neutrino daughter, not primary lepton",400,0,40);
 TH1D * hothernudauP_charged = new TH1D("hothernudauP_charged","Charged nu product, not primary lepton, not charm",400,0,40);
 TH1D * hothernudauP_charm = new TH1D("hothernudauP_charm","Charmed hadron, nu product",400,0,40);
 TH1D * hothernudauPdgCode = new TH1D("hothernudauPdgCode","PdgCode of neutrino daughters, apart the primary lepton",6000,-3000.,3000.);

 TH1D * hmuP = new TH1D("hmuP","Muon momentum;P[GeV/c]",400,0,400);
 TH1D * hmuTheta = new TH1D("hmuTheta","Theta angle of muonsw;#theta[rad]",350,0.,3.5);
 TH1D * hmuPhi = new TH1D("hmuPhi","Phi angle of muonsw;#phi[rad]",700,-3.5,3.5);

 TH1D *h_tot_charged_hadron_P = new TH1D("h_tot_charged_hadron_P","total hadron momentum",400,0,400);

 //2D Theta_P
 TH2D * hmuTheta_P = new TH2D("hmuTheta_P","Theta angle of muon vs momentum;P[GeV/c];#theta[rad]",400,0,400,350,0.,3.5);

 //mu- histograms
 TH2D *hdy_dz = new TH2D("hdy_dz","At the end of spectrometer, Y distance from vy with respect to distance from vz;dz[cm];dy[cm]",30,300.,600.,60,-300.,300.);
 //TH1D *hdeltaTX_muminus = new TH1D("hdeltaTX_muminus","TX difference;DeltaTX",400,-0.02,0.02);
 TH1D *hdeltaTY_muminus = new TH1D("hdeltaTY_muminus","y projection sagitta negative muons;SY[cm]",400,-20.,20.);
 //mu+ histograms
 //TH1D *hdeltaTX_muplus = new TH1D("hdeltaTX_muplus","TX difference;DeltaTX",400,-0.02,0.02);
 TH1D *hdeltaTY_muplus = new TH1D("hdeltaTY_muplus","y projection of sagitta positive muons;SY[cm]",400,-20.,20.);

 TH1D *hdeltaTX = new TH1D("hdeltaTX","x projection of sagitta;SX[cm]",40,-0.2,0.2);

 TProfile *prof_deltaTY_pzy = new TProfile("prof_deltaTY_pzy","momentum vs angular difference;1/dTy;pzy[GeV/c];",200,0,200.,0,400); //deltaTY here is considered as abs
 TProfile *prof_deltaTY_p = new TProfile("prof_deltaTY_p","momentum vs angular difference;1/dTy;p[GeV/c]",200,0,200.,0,400); //deltaTY here is considered as abs

 TF1 *fang_pzy = new TF1("fang_pzy","pol1",0,200.);

 TH1D *hres_pzy = new TH1D("hres_pzy","resolution of transverse momentum;res_pzy{GeV/c}",200,-1.,1.);
 //fang_pzy->SetParameter(0,0.);
 //fang_pzy->SetParameter(0,1.4576);
 //fang_pzy->SetParameter(1,1.8754);

 fang_pzy->SetParameter(0,1.47952);
 fang_pzy->SetParameter(1,1.80416);

 const int nstations = 4;
 TH2D *hxy[nstations];
 TH2D *hxy_acceptance[nstations];
 //initialize histos
 for (int istation = 0; istation < nstations; istation++){
  hxy[istation] = new TH2D(Form("hxy[%i]",istation),Form("Muon distribution in station %i;x[cm];y[cm]",istation),400,-200,200,400,-200,200);
  hxy_acceptance[istation] = new TH2D(Form("hxy_iacceptance[%i]",istation),Form("Muon distribution in station %i acceptance studu;x[cm];y[cm]",istation),
    (int)(2*dx_acceptance),-dx_acceptance,dx_acceptance,(int)(2*dy_acceptance),-dy_acceptance,dy_acceptance);
 }
 Double_t XSpectrometer[nstations], YSpectrometer[nstations], ZSpectrometer[nstations];
 Bool_t InSpectrometer[nstations]; //check if the muon passed the station

 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;

 int totalweight = 0; //for integral normalization

 const double init_xspectro = -99999.; //big value for check if they are present

 const float sigma_delta_theta = 6.91005e-04; //as measured with previous iteration

 double chargeeff = 0.;

 double tot_charged_hadron_P;

 for(int ientry = 0;ientry<nentries;ientry++){
  tot_charged_hadron_P = 0.; //resetting sum
  //only even number
  if (ientry%2==eventesclusion) continue;
  //array initialization
  for (int istation = 0; istation < nstations; istation++){
    XSpectrometer[istation] = init_xspectro;
    YSpectrometer[istation] = init_xspectro;
    ZSpectrometer[istation] = init_xspectro;
    InSpectrometer[istation] = false;
  }
  reader.SetEntry(ientry);
  Double_t weight = tracks[0].GetWeight();
  totalweight += weight;
  Double_t vz = tracks[0].GetStartZ();
  Double_t vy = tracks[0].GetStartY();
  if (tracks.GetSize() < 2){
    cout<<"WARNING: Only neutrino in event "<<ientry<<endl;
    continue;
  }

  //Double_t weight = tracks[1].GetWeight();
  //loop on tracks
  int itrack = 0;
  for (const ShipMCTrack& track: tracks){
    bool ischarm = false;
    int pdgcode = track.GetPdgCode();
    //check if the track is a charmed hadron (pdgcode belongs to list)
    auto searchresult = std::find(begin(charmpdglist),end(charmpdglist),TMath::Abs(pdgcode));
    if (searchresult != end(charmpdglist)) ischarm = true;
    //end of the check
    double charge = GetParticleCharge(pdgcode, pdg);
    if (itrack == 1){ 
    //only looking at the initial muon
     if (TMath::Abs(pdgcode)==13){ //saving initial muon information
     profmuE_nuE->Fill(track.GetEnergy(), track.GetEnergy()/tracks[0].GetEnergy()); //no weight here, it is a tprofile, so we want the correlation
     hmuP->Fill(track.GetP(),weight);

     TVector3 *muonmom = new TVector3(track.GetPx(),track.GetPy(),track.GetPz());
     hmuTheta->Fill(muonmom->Theta(),weight);
     hmuPhi->Fill(muonmom->Phi(),weight);
    //2D plot theta at different momenta
     hmuTheta_P->Fill(track.GetP(),muonmom->Theta(),weight);
    //Double_t mytheta = TMath::ACos(tracks[1].GetPz()/tracks[1].GetP());
    //cout<<mytheta<<" "<<muonmom->Theta()<<endl; //Naturally, they are the same. Keep here just as a Memento of my paranoia
    //storing muon information
   }
    else{
    cout<<"WARNING: Track 1 is not muon neutrino in event "<<ientry<<endl;
    continue; //muon missing in this event (too low kin energy)
   } 
  }
    else{
    //is it a charged neutrino daughter (not the primary lepton)?
     if (track.GetMotherId()==0){
      hothernudauP->Fill(track.GetP(),weight);
      if (TMath::Abs(charge) > 0 && !ischarm) hothernudauP_charged->Fill(track.GetP(),weight);
      if (ischarm) hothernudauP_charm->Fill(track.GetP(),weight);
      if (TMath::Abs(charge) > 0 ) tot_charged_hadron_P += track.GetP(); //we want to plot the total hadron energy;
      hothernudauPdgCode->Fill(track.GetPdgCode(),weight);
     }
   }
   itrack++; //trackID counter 
  }//end loop over tracks
  h_tot_charged_hadron_P->Fill(tot_charged_hadron_P, weight);
 //first, loop over all hits to find positions of muon in spectrometer
 for (const ShipRpcPoint& rpcpoint: rpcpoints){
  int nstation = rpcpoint.GetDetectorID() - 1; //from 0 to 3 
  if (rpcpoint.GetTrackID()==1 && !InSpectrometer[nstation]){ //we fill the arrays with the muon neutrino positions
  
   XSpectrometer[nstation] = rpcpoint.GetX();
   YSpectrometer[nstation] = rpcpoint.GetY();
   ZSpectrometer[nstation] = rpcpoint.GetZ();
   InSpectrometer[nstation] = true;

   hxy[nstation]->Fill(XSpectrometer[nstation],YSpectrometer[nstation], weight);
   hxy_acceptance[nstation]->Fill(XSpectrometer[nstation],YSpectrometer[nstation], weight);
   //applying smearing to X and Y
   XSpectrometer[nstation] = XSpectrometer[nstation] + randomgen->Gaus(0, spectro_posres);
   YSpectrometer[nstation] = YSpectrometer[nstation] + randomgen->Gaus(0, spectro_posres);
  }//end muon track condition
 }//end rpc hit loop 

  //check if hits in all stations are present
  if (InSpectrometer[0]&&InSpectrometer[1]&&InSpectrometer[2]&&InSpectrometer[3]){
    //auto DeltaT = BendingSpectrometer(XSpectrometer, YSpectrometer, ZSpectrometer);
    auto DeltaT = SagittaSpectrometer(XSpectrometer, YSpectrometer, ZSpectrometer);
    hdeltaTX->Fill(DeltaT[0],weight);
    if(tracks[1].GetPdgCode()==13){ //negative muon
    hdeltaTY_muminus->Fill(DeltaT[1],weight);
    }
    else if(tracks[1].GetPdgCode()==-13) { //positive muon
    hdeltaTY_muplus->Fill(DeltaT[1],weight);
    }

    if (TMath::Abs(DeltaT[1]) > 3 * sigma_delta_theta) chargeeff+=weight; //as usual, requiring 3 sigma for charge separation

    Double_t pzy = TMath::Sqrt(pow(tracks[1].GetPy(),2) + pow(tracks[1].GetPz(),2));//opposite to magnetic field
    prof_deltaTY_p->Fill(1./TMath::Abs(DeltaT[1]),tracks[1].GetP());
    prof_deltaTY_pzy->Fill(1./TMath::Abs(DeltaT[1]),pzy);

    Double_t meas_pzy = fang_pzy->Eval(1./TMath::Abs(DeltaT[1]));
    Double_t res_pzy = (meas_pzy - pzy)/pzy;

    hres_pzy->Fill(res_pzy,weight);

    hdy_dz->Fill(ZSpectrometer[3]-vz,YSpectrometer[3]-vy,weight); //How far it goes from vertex position?
  
  }
  
 } 

 //**CHECKING ACCEPTANCE**//
 cout<<"Fraction of hits arriving at station 2: "<<hxy[1]->Integral()/totalweight<<endl;
 cout<<"Subfraction of hits in acceptance station 2: "<<hxy_acceptance[1]->Integral()/hxy[1]->Integral()<<endl;
 cout<<"Total acceptance station 2: "<<hxy_acceptance[1]->Integral()/totalweight<<endl;

 cout<<"Fraction of hits arriving at station 4: "<<hxy[3]->Integral()/totalweight<<endl;
 cout<<"Subfraction of hits in acceptance: "<<hxy_acceptance[3]->Integral()/hxy[3]->Integral()<<endl;
 cout<<"Total acceptance: "<<hxy_acceptance[3]->Integral()/totalweight<<endl;

 //**CHECKING CHARGE IDENTIFICATION**//
 cout<<"Fraction of weighted muon events with charge identified: "<<chargeeff/hxy[3]->Integral()<<endl;
 //plotting normalized distributions
 TCanvas *c = new TCanvas();
 hmuP->Scale(1./hmuP->Integral());
 hmuP->Draw("histo");
 //comparing lepton and total (charged only) hadron momentum
 h_tot_charged_hadron_P->Scale(1./h_tot_charged_hadron_P->Integral());
 h_tot_charged_hadron_P->SetLineColor(kRed);
 h_tot_charged_hadron_P->Draw("histo&&SAME");
 c->BuildLegend();

 TCanvas *cangles = new TCanvas();
 hmuTheta->Scale(1./hmuTheta->Integral());
 hmuTheta->Draw();

 TCanvas *cprof_muEnuE = new TCanvas();
 profmuE_nuE->Draw();

 TCanvas *cothernudau = new TCanvas();
 hothernudauP->Scale(1./hothernudauP->Integral());
 //hothernudauP->Draw();
 hothernudauP_charged->Scale(1./hothernudauP_charged->Integral());
 hothernudauP_charged->SetLineColor(kRed);
 //hothernudauP_charm->SetLineColor(kYellow);
 hothernudauP_charged->Draw();
 hothernudauP_charm->Scale(1./hothernudauP_charm->Integral());
 hothernudauP_charm->Draw("SAME");
 cothernudau->BuildLegend();
 
 TCanvas *cotherpdgcode = new TCanvas();
 hothernudauPdgCode->Draw();

 //TCanvas *cdy_vz = new TCanvas("cdy_vz","Distance as a function of vertex position");
 //hdy_dz->Scale(1./hdy_dz->Integral());
 //hdy_dz->Draw("COLZ");

 TCanvas *cxy_up = new TCanvas("cxy_up","upstream distribution of muons",800,800);
 hxy[0]->Scale(1./hxy[0]->Integral());
 hxy[0]->Draw();

 TCanvas *cxy_middle = new TCanvas("cxy_down","downstream distribution of muons",800,800);
 hxy[1]->Scale(1./hxy[1]->Integral());
 hxy[1]->Draw();

 TCanvas *cxy_down = new TCanvas("cxy_down","downstream distribution of muons",800,800);
 hxy[3]->Scale(1./hxy[3]->Integral());
 hxy[3]->Draw();

 TCanvas *cspectroangles = new TCanvas("cspectroangles","Angular differences from spectrometer");
 cspectroangles->Divide(2,1);
 cspectroangles->cd(1);
 hdeltaTX->Scale(1./hdeltaTX->Integral());
 hdeltaTX->Fit("gaus");
 hdeltaTX->Draw("histo");
 hdeltaTX->GetFunction("gaus")->Draw("SAME");
 cspectroangles->cd(2);
 hdeltaTY_muminus->Scale(1./hdeltaTY_muminus->Integral());
 hdeltaTY_muminus->SetLineColor(kBlue);
 hdeltaTY_muminus->Draw("histo");
 //hdeltaTY_muplus->Scale(1./hdeltaTY_muplus->Integral());
 //hdeltaTY_muplus->SetLineColor(kRed);
 //hdeltaTY_muplus->Draw("histo && SAMES");
 cspectroangles->GetPad(2)->BuildLegend();
 hdeltaTY_muminus->SetTitle("");

 //TCanvas *cprof_p = new TCanvas();
 //prof_deltaTY_p->Draw();

 TCanvas *cprof_pT = new TCanvas();
 prof_deltaTY_pzy->SetMarkerStyle(kFullCircle);
 prof_deltaTY_pzy->Draw();
 //fang_pzy->FixParameter(0,0);
 prof_deltaTY_pzy->Fit(fang_pzy);

 TCanvas *chres_pzy = new TCanvas();
 hres_pzy->Scale(1./hres_pzy->Integral());
 hres_pzy->Draw("histo");
 hres_pzy->Fit("gaus");
 hres_pzy->GetFunction("gaus")->Draw("SAME");

 TCanvas *cmutheta_P = new TCanvas();
 hmuTheta_P->Scale(1./hmuTheta_P->Integral());
 hmuTheta_P->Draw("COLZ");

 TCanvas *ctothadronP = new TCanvas();
 h_tot_charged_hadron_P->Scale(1./h_tot_charged_hadron_P->Integral());
 h_tot_charged_hadron_P->Draw();
}
