//we need to know how many showers are produced per event in our brick and their energy/molteplicity
void shower_loop(TString inputfile){

 TFile *file = TFile::Open(inputfile.Data()); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 //defining histograms   
 TH1D * hstartz = new TH1D("hstartz","Startz position;z[cm]",80,121.6,125.6);
 TH1D * henergy = new TH1D("henergy","Distribution in energy;E[GeV]",100,0,50);
 TH1I * hpdgcode = new TH1I("hpdgcode","pdgcode",50,-25,25);
 TH1I * hnshowers = new TH1I("hnshowers","Number of shower per event;nshowers",200,0,200);

 double endtarget = 125.6;

 const int nentries = 10000;
// const int nentries = 10;

 int motherid, procid;
 map<int, bool> shower_initiators;

 

 int pairprocid = 5;
 int bremprocid = 8;

 cout<<"Starting loop over "<<nentries<<endl;

 TDatabasePDG *pdg = TDatabasePDG::Instance();
 bool fromshower;
 float showerstartz;
 for (int ientry = 0; ientry < nentries; ientry++){
  if (ientry % 1000==0) cout<<"Arrived at entry "<<ientry<<endl;
  reader.SetEntry(ientry);
  //resetting container of showers
  shower_initiators.clear();
  //starting loop over tracks
  for (int itrk = 0; itrk < tracks.GetSize(); itrk++){
   fromshower = false;
   int trackID = itrk;
   ShipMCTrack track = tracks[itrk];
   if (track.GetStartZ() < endtarget){ //we are interested only in showers in our target
    while (track.GetProcID()==pairprocid || track.GetProcID()==bremprocid){ //track from shower, go to mothers until the origin of the shower
     fromshower = true;
     trackID = track.GetMotherId();
     showerstartz = track.GetStartZ();
     track = tracks[trackID]; //passing to the mother
    }
   }//end condition about z position
   if (fromshower && (shower_initiators[trackID] == false)){ //filling histograms
    shower_initiators[trackID] = true; //found
    henergy->Fill(track.GetEnergy());
    hpdgcode->Fill(track.GetPdgCode());
    hstartz->Fill(showerstartz);
   } //end condition of found shower (it requires also the z position to be true)
  }//end loop over tracks
   hnshowers->Fill(shower_initiators.size());
 }//end loop over events
 TCanvas *cenergy = new TCanvas();
 henergy->Draw();
 TCanvas * cstartz = new TCanvas();
 hstartz->Draw();
 TCanvas * cnshowers = new TCanvas();
 hnshowers->Draw();
 TCanvas * cpdgcode = new TCanvas();
 hpdgcode->Draw();
}
