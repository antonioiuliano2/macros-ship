//fast drawing, no need to loop over all tracks in file. Readout is done manually. Does not work for more than one track.
void drawEDAtrack(int trackID = 129){

 EdbPVRec *ali = new EdbPVRec();
 TObjArray *drawntracks = new TObjArray(100);
  //reading track file and setting branches
 TFile *trackfile = TFile::Open("linked_tracks.root");
 TTree *tracks = (TTree*)trackfile->Get("tracks");

  Int_t   trid=0;
  Int_t   nseg=0;
  Int_t   npl=0;
  Int_t   n0=0;
  Float_t xv=0.;
  Float_t yv=0.;

  TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
  TClonesArray *segf = new TClonesArray("EdbSegP", 60);
  EdbSegP *trk=0;
  EdbSegP *s1=0;

  tracks->SetBranchAddress("trid", &trid);
  tracks->SetBranchAddress("nseg", &nseg);
  tracks->SetBranchAddress("npl", &npl);
  tracks->SetBranchAddress("n0", &n0);
  tracks->SetBranchAddress("xv", &xv);
  tracks->SetBranchAddress("yv", &yv);
  tracks->SetBranchAddress("sf", &segf);
  tracks->SetBranchAddress("s",  &seg);
  tracks->SetBranchAddress("t.", &trk);

  //Getting Entry

  tracks->BuildIndex("trid");
  tracks->GetEntryWithIndex(trackID);

  EdbTrackP *tr1 = new EdbTrackP();
  ((EdbSegP*)tr1)->Copy(*trk);
  tr1->SetM(0.139);                 //TODO

  for(int i=0; i<nseg; i++) {
      s1 = (EdbSegP*)(seg->At(i));

      tr1->AddSegment( s1 );
      tr1->AddSegmentF( new EdbSegP(*((EdbSegP*)(segf->At(i)))) );
    }
    tr1->SetSegmentsTrack(tr1->ID());
    tr1->SetCounters();
    //tr1->FitTrackKFS(true);

 //adding track to drawing set and building EDA
 if(tr1->ID()==trackID) drawntracks->Add(tr1);
 ali->eTracks = drawntracks;
 //ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->Run();
}
