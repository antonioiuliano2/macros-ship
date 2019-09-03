void drawEDAvertex(int vID = 0){
 const int nvertices = 1;
 int vertexlist[nvertices] = {2947};
 int trackID = 166641;

 TFile * inputfile = TFile::Open("vertices_firstquarter.root");
 EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");


 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);


 for (int ivtx = 0; ivtx < nvertices; ivtx++){
  int vID = vertexlist[ivtx];
  EdbVertex *vertex = (EdbVertex*) vertexrec->eVTX->At(vID);

  drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
  for (int itrk = 0; itrk < vertex->N(); itrk++){
     EdbTrackP* track =  vertex->GetTrack(itrk);
     drawntracks->Add(track);
  }
 }
 EdbPVRec *ali = new EdbPVRec();

  //reading track file and setting branches
  TFile *trackfile = TFile::Open("verticesandtracks.root");
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
 cout<<"Test: "<<tr1->N()<<endl;
 drawntracks->Add(tr1);
 ali->eTracks = drawntracks;
 ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->Run();
}
