//loop over stand-alone G4 pot background (27 Gennaio)
using namespace ROOT;
void G4_bkgloop(TString inputfile){
 TFile *file = TFile::Open(inputfile.Data()); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 TFile * outputfile = new TFile("histos_primary_bkg.root","RECREATE");

 //defining histograms   
 TH1D * hstartz = new TH1D("hstartz","Startz position;z[#mum]",100,120,130);
 TH1D * henergy = new TH1D("henergy","Distribution in energy;E[GeV]",400,0,400);
 TH1I * hnprimaries = new TH1I("hnprimaries", "Number of tracks at primary vertex",100,0,100);

 double endtarget = 125.6;

 const int nentries = reader.GetEntries();

 int motherid, procid;
 int nprimaries;
 RVec<double> energy_primaries;
 double mindeltaE = 10.;
 double startz;

 cout<<"Starting loop over "<<nentries<<endl;

 TDatabasePDG *pdg = TDatabasePDG::Instance();
 for (int ientry = 0; ientry < nentries; ientry++){
  if (ientry % 10000==0) cout<<"Arrived at entry "<<ientry<<endl;
  reader.SetEntry(ientry);
  //starting loop over tracks
  nprimaries = 0;
  startz = 0.;
  energy_primaries.clear();
  int pdgcode;
  double primaryenergy = tracks[0].GetEnergy();
  bool interactingintarget = false;
  for (const ShipMCTrack &MCTrack:tracks){

   motherid = MCTrack.GetMotherId();
   procid = MCTrack.GetProcID();

   if (procid == 23 && motherid==0){ // primary daughter selection
     startz = MCTrack.GetStartZ();
     pdgcode = MCTrack.GetPdgCode();
     int charge = 0;
     if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
     if(startz<endtarget&&TMath::Abs(charge)>0){ 
      nprimaries++;
      energy_primaries.push_back(MCTrack.GetEnergy());    
     }
    }
  } //end loop over tracks
  if (nprimaries > 0){
    if (primaryenergy - Max(energy_primaries)>mindeltaE){
     hnprimaries->Fill(nprimaries);
     hstartz->Fill(startz);
     for (const double energy: energy_primaries) henergy->Fill(energy);
    }
    //if (nprimaries < 3) cout<<"Suspect event "<<ientry<<endl;
   }
 }//end loop over events
 outputfile->cd();
 hstartz->Write();
 henergy->Write();
 hnprimaries->Write();
 outputfile->Close();
}
