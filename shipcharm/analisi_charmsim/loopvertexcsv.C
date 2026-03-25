void loopvertexcsv(){

  TTree *MCtree = new TTree("MCtree","MC categories after FEDRA reco");
  MCtree->ReadFile("MC_vertexlist_withsmearing.csv","ntracks/I:ivtx/I:itrk/I:MCEventID/I:MCTrackID/I:MCMotherID/I:predmolt/I:preddecaylength/F:quantity/I:vx/F:vy/F:vz/F:topology/I");

  //setting the reader for the tree
  TTreeReader treereader = TTreeReader(MCtree);
  TTreeReaderValue<Int_t> ntracks(treereader,"ntracks");
  TTreeReaderValue<Int_t> predmolt(treereader,"predmolt");
  TTreeReaderValue<Int_t> ivtx(treereader,"ivtx");
  TTreeReaderValue<Int_t> itrk(treereader,"itrk");
  TTreeReaderValue<Int_t> MCEventID(treereader,"MCEventID");
  TTreeReaderValue<Int_t> MCTrackID(treereader,"MCTrackID");
  TTreeReaderValue<Int_t> MCMotherID(treereader,"MCMotherID");
  TTreeReaderValue<Int_t> topology(treereader,"topology");
  
  TTreeReaderValue<Float_t> preddecaylength(treereader,"preddecaylength");
  const int nentries = treereader.GetEntries();
  //simulation tree
  Double_t charmenergy;
  TFile *simfile = TFile::Open("ship.conical.Pythia8CharmOnly-TGeant4.root");
  TTreeReader simreader("cbmsim",simfile);
  TTreeReaderArray<ShipMCTrack> tracks (simreader,"MCTrack");
 
  for (int ientry = 0; ientry<nentries;ientry++){
    treereader.Next();
    cout<<*topology<<" "<<*itrk<<endl;
    if (*topology == 1) charmenergy = 0.; //not charm daughter, but primary track
    else{
      simreader.SetEntry(*MCEventID);
      ShipMCTrack charmtrack = tracks[*MCMotherID];
      charmenergy = charmtrack.GetEnergy();
    }
  }
}
