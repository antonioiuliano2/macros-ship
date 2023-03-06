//reading numu and numubar simulation, AFTER they have cut within target acceptance with cutsim.py
//verify muon charge separation from angular smearing along y
void muonspectrometer_testsimulation(){

 const int eventesclusion = 2; //0 only odd, 1 only even, 2 all;
 const float dx_acceptance = 60.; //how many we lose by reducing our station size?
 const float dy_acceptance = 60.; 

 TString prefix("root:://eosuser.cern.ch/");//for ROOTXD
 TString simpath_mu("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/numu_CCDIS_2023_03_04_targetmovedupstream_spectromag/");
 TString simpath_mu_bar("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/numu_bar_CCDIS_2023_03_04_targetmovedupstream_spectromag/");

 TRandom3 *randomgen  = new TRandom3(); //0; seed changes everytime, no seed, default to 4357 

 const double spectro_posres = 100.*1e-4; //100 micron, for smearing 
 
 TChain *simchain = new TChain("cbmsim");
 simchain->Add((prefix+simpath_mu+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 simchain->Add((prefix+simpath_mu_bar+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 
 TTreeReader reader(simchain);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 TH1D * hmuE = new TH1D("hmuE","Muon energy;E[GeV]",400,0,400);
 //mu- histograms
 TH2D *hdy_dz = new TH2D("hdy_dz","At the end of spectrometer, Y distance from vy with respect to distance from vz;dz[cm];dy[cm]",30,300.,600.,60,-300.,300.);
 //TH1D *hdeltaTX_muminus = new TH1D("hdeltaTX_muminus","TX difference;DeltaTX",400,-0.02,0.02);
 TH1D *hdeltaTY_muminus = new TH1D("hdeltaTY_muminus","TY difference;DeltaTY",4000,-2.,2.);
 //mu+ histograms
 //TH1D *hdeltaTX_muplus = new TH1D("hdeltaTX_muplus","TX difference;DeltaTX",400,-0.02,0.02);
 TH1D *hdeltaTY_muplus = new TH1D("hdeltaTY_muplus","TY difference;DeltaTY",4000,-2.,2.);

 TH1D *hdeltaTX = new TH1D("hdeltaTX","TX difference;DeltaTX",400,-0.02,0.02);

 TProfile *prof_deltaTY_pzy = new TProfile("prof_deltaTY_pzy","momentum vs angular difference;pzy[GeV/c];1/dTy",200,0,200.,0,400); //deltaTY here is considered as abs
 TProfile *prof_deltaTY_p = new TProfile("prof_deltaTY_p","momentum vs angular difference;p[GeV/c];1/dTy",200,0,200.,0,400); //deltaTY here is considered as abs

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


 for(int ientry = 0;ientry<nentries;ientry++){
  
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
  if (TMath::Abs(tracks[1].GetPdgCode())==13){ //saving initial muon information

    hmuE->Fill(tracks[1].GetEnergy(),weight);
    //storing muon information
   }
   else{
    cout<<"WARNING: Track 1 is not muon neutrino in event "<<ientry<<endl;
    continue; //muon missing in this event (too low kin energy)
   } 
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
    //0 and 1: upstream stations;
    Double_t TXup = (XSpectrometer[1] - XSpectrometer[0])/(ZSpectrometer[1] - ZSpectrometer[0]);
    Double_t TYup = (YSpectrometer[1] - YSpectrometer[0])/(ZSpectrometer[1] - ZSpectrometer[0]);
    //2 and 3: downstream stations;
    Double_t TXdown = (XSpectrometer[3] - XSpectrometer[2])/(ZSpectrometer[3] - ZSpectrometer[2]);
    Double_t TYdown = (YSpectrometer[3] - YSpectrometer[2])/(ZSpectrometer[3] - ZSpectrometer[2]);

    Double_t DeltaTX = TXdown - TXup;
    Double_t DeltaTY = TYdown - TYup;
    hdeltaTX->Fill(DeltaTX,weight);
    if(tracks[1].GetPdgCode()==13){
    //hdeltaTX_muminus->Fill(DeltaTX,weight);
    hdeltaTY_muminus->Fill(DeltaTY,weight);
    }
    else{
    //hdeltaTX_muplus->Fill(DeltaTX,weight);
    hdeltaTY_muplus->Fill(DeltaTY,weight);
    }
    Double_t pzy = TMath::Sqrt(pow(tracks[1].GetPy(),2) + pow(tracks[1].GetPz(),2));//opposite to magnetic field
    prof_deltaTY_p->Fill(1./TMath::Abs(DeltaTY),tracks[1].GetP());
    prof_deltaTY_pzy->Fill(1./TMath::Abs(DeltaTY),pzy);

    Double_t meas_pzy = fang_pzy->Eval(1./TMath::Abs(DeltaTY));
    Double_t res_pzy = (meas_pzy - pzy)/pzy;

    hres_pzy->Fill(res_pzy,weight);

    hdy_dz->Fill(ZSpectrometer[3]-vz,YSpectrometer[3]-vy,weight); //How far it goes from vertex position?
  
  }
  
 } 

 //**CHECKING ACCEPTANCE**//
 cout<<"Fraction of hits arriving at station 4: "<<hxy[3]->Integral()/totalweight<<endl;
 cout<<"Subfraction of hits in acceptance: "<<hxy_acceptance[3]->Integral()/hxy[3]->Integral()<<endl;
 cout<<"Total acceptance: "<<hxy_acceptance[3]->Integral()/totalweight<<endl;

 //plotting normalized distributions
 TCanvas *c = new TCanvas();
 hmuE->Scale(1./hmuE->Integral());
 hmuE->Draw("histo");

 TCanvas *cdy_vz = new TCanvas("cdy_vz","Distance as a function of vertex position");
 hdy_dz->Scale(1./hdy_dz->Integral());
 hdy_dz->Draw("COLZ");

 TCanvas *cxy_up = new TCanvas("cxy_up","upstream distribution of muons",800,800);
 hxy[0]->Scale(1./hxy[0]->Integral());
 hxy[0]->Draw();

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
 hdeltaTY_muplus->Scale(1./hdeltaTY_muplus->Integral());
 hdeltaTY_muplus->SetLineColor(kRed);
 hdeltaTY_muplus->Draw("histo && SAMES");

 TCanvas *cprof_pT = new TCanvas();
 prof_deltaTY_p->Draw();

 TCanvas *cprof_p = new TCanvas();
 prof_deltaTY_pzy->SetMarkerStyle(kFullCircle);
 prof_deltaTY_pzy->Draw();
 //fang_pzy->FixParameter(0,0);
 prof_deltaTY_pzy->Fit(fang_pzy);

 TCanvas *chres_pzy = new TCanvas();
 hres_pzy->Scale(1./hres_pzy->Integral());
 hres_pzy->Draw("histo");
 hres_pzy->Fit("gaus");
 hres_pzy->GetFunction("gaus")->Draw("SAME");
}
