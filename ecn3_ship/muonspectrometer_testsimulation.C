//reading numu simulation, AFTER they have cut within target acceptance
void muonspectrometer_testsimulation(){
 
 TFile *inputfile = TFile::Open("inECC_ship.conical.Genie-TGeant4.root");
 TTree *simtree = (TTree*) inputfile->Get("cbmsim");
 
 TTreeReader reader(simtree);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 TH1D * hmuE = new TH1D("hmuE","Muon energy;E[GeV]",400,0,400);

 TH2D *hxy_upstream = new TH2D("hxy_upstream","Muon distribution upstream;x[cm];y[cm]",1000,-500,500,1000,-500,500);
 TH2D *hxy_downstream = new TH2D("hxy_downstream","Muon distribution downstream;x[cm];y[cm]",1000,-500,500,1000,-500,500);

 const int nstations = 4;
 Double_t XSpectrometer[nstations], YSpectrometer[nstations], ZSpectrometer[nstations];

 const double 

 

 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;

 for(int ientry = 0;ientry<nentries;ientry++){   
  reader.SetEntry(ientry);
  Double_t weight = tracks[1].GetWeight();
  if (tracks[1].GetPdgCode()==13){ //saving initial muon information

    hmuE->Fill(tracks[1].GetEnergy(),weight);
    //storing muon information
   } 

 for (const ShipRpcPoint& rpcpoint: rpcpoints){
  if (rpcpoint.GetTrackID()==1){ //we only study the muon here
   if (rpcpoint.GetDetectorID()==1){

    hxy_upstream->Fill(rpcpoint.GetX(), rpcpoint.GetY(), weight);

   }

  if (rpcpoint.GetDetectorID()==4){

    hxy_downstream->Fill(rpcpoint.GetX(), rpcpoint.GetY(), weight);

   }
  }//end muon track condition
 }//end rpc hit loop 
  
  
 } 
 TCanvas *c = new TCanvas();
 hmuE->Scale(1./hmuE->Integral());
 hmuE->Draw("histo");

 TCanvas *cxy_up = new TCanvas("cxy_up","upstream distribution of muons",800,800);
 hxy_upstream->Draw();

 TCanvas *cxy_down = new TCanvas("cxy_down","downstream distribution of muons",800,800);
 hxy_downstream->Draw();

 


}
