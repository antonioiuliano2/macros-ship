//reading numu simulation, AFTER they have cut within target acceptance
void muonspectrometer_testsimulation(){

 
 TString simpath("root:://eosuser.cern.ch//eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/numuCCDIS_dummyspectrometer_2023_February_27/magspectro_performance/");

 TRandom3 *randomgen  = new TRandom3(); //0; seed changes everytime, no seed, default to 4357 

 const double spectro_posres = 100.*1e-4; //100 micron, for smearing 
 
 TFile *inputfile = TFile::Open((simpath+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 TTree *simtree = (TTree*) inputfile->Get("cbmsim");
 
 TTreeReader reader(simtree);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 TH1D * hmuE = new TH1D("hmuE","Muon energy;E[GeV]",400,0,400);

 TH1D *hdeltaTX = new TH1D("hdeltaTX","TX difference;DeltaTX",400,-0.02,0.02);
 TH1D *hdeltaTY = new TH1D("hdeltaTY","TY difference;DeltaTY",4000,-2.,2.);

 const int nstations = 4;
 TH2D *hxy[nstations];
 //initialize histos
 for (int istation = 0; istation < nstations; istation++){
  hxy[istation] = new TH2D(Form("hxy[%i]",istation),Form("Muon distribution in station %i;x[cm];y[cm]",istation),1000,-500,500,1000,-500,500);
 }
 Double_t XSpectrometer[nstations], YSpectrometer[nstations], ZSpectrometer[nstations];
 Bool_t InSpectrometer[nstations]; //check if the muon passed the station

 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;

 const double init_xspectro = -99999.; //big value for check if they are present


 for(int ientry = 0;ientry<nentries;ientry++){
  //array initialization
  for (int istation = 0; istation < nstations; istation++){
    XSpectrometer[istation] = init_xspectro;
    YSpectrometer[istation] = init_xspectro;
    ZSpectrometer[istation] = init_xspectro;
    InSpectrometer[istation] = false;
  }
  reader.SetEntry(ientry);
  if (tracks.GetSize() < 2){
    cout<<"WARNING: Only neutrino in event "<<ientry<<endl;
    continue;
  }

  Double_t weight = tracks[1].GetWeight();
  if (tracks[1].GetPdgCode()==13){ //saving initial muon information

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
    hdeltaTY->Fill(DeltaTY,weight);
  }
  
 } 

 //plotting
 TCanvas *c = new TCanvas();
 hmuE->Scale(1./hmuE->Integral());
 hmuE->Draw("histo");

 TCanvas *cxy_up = new TCanvas("cxy_up","upstream distribution of muons",800,800);
 hxy[0]->Draw();

 TCanvas *cxy_down = new TCanvas("cxy_down","downstream distribution of muons",800,800);
 hxy[3]->Draw();

 TCanvas *cspectroangles = new TCanvas("cspectroangles","Angular differences from spectrometer");
 cspectroangles->Divide(1,2);
 cspectroangles->cd(1);
 hdeltaTX->Fit("gaus");
 hdeltaTX->Draw("histo");
 hdeltaTX->GetFunction("gaus")->Draw("SAME");
 cspectroangles->cd(2);
 hdeltaTY->Draw("histo");


}
