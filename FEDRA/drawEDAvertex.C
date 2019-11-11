//display for vertices, patterns are not drawn

EdbVertex* GetVertexFromTree( EdbPVRec &ali, const char     *fname, const int vertexID );

void drawEDAvertex(bool newversion = true, TString vertexfilename= "vertextree.root"){
 const int nvertices = 1;
 int vertexlist[nvertices] = {139653, 40204};
 int vertexcolors[nvertices] = {kRed,kGreen};

 TFile * inputfile = TFile::Open(vertexfilename.Data());
 EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");


 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);

// EdbTrackP *specialtrack = new EdbTrackP();

 EdbPVRec *ali = new EdbPVRec();

 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = vertexlist[i];
  
  EdbVertex *vertex = 0;

  if (newversion) vertex = GetVertexFromTree(*ali,vertexfilename,vID);
  else vertex = (EdbVertex*) vertexrec->eVTX->At(vID);


  drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
  for (int itrk = 0; itrk < vertex->N(); itrk++){
     EdbTrackP* track =  vertex->GetTrack(itrk);     
     for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[i]); // to color them differently
 //    for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[i]); // to color them differently
     drawntracks->Add(track);
//     specialtrack = track;
  }
 }

 ali->eTracks = drawntracks;
 ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
/* eda->GetTrackSet("TS")->RemoveTrack(specialtrack); //charm colored differently
 eda->GetTrackSet("SB")->SetTrackAttribute(4);
 eda->GetTrackSet("SB")->AddTrack(specialtrack);*/
 eda->Run();
}

EdbVertex* GetVertexFromTree( EdbPVRec &ali, const char     *fname, const int vertexID )
{
  TFile f(fname);
 // if(f.IsZombie()) { Log(1,"EdbDataProc::ReadVertexTree","Error open file %s", fname);  return 0; }
  
  TTree *vtx = (TTree*)f.Get("vtx");

  //vertex variables
  Float_t vx, vy, vz;
  Int_t vID = 0;
  Int_t n = 0;
  Int_t flag = 0;
  //track variables

  const Int_t maxdim = 1000; //maximum number of tracks associated to a vertex that can be saved in the tree
  //big arrays for containers of track variables
  Int_t TrackID[maxdim];
  Int_t nholes[maxdim];
  Int_t maxgap[maxdim];
  Int_t nseg[maxdim];
  Float_t TX[maxdim];
  Float_t TY[maxdim];
  Float_t impactparameter[maxdim];
  Int_t incoming[maxdim]; //coming from the vertex
  //MC information
  Int_t MCEventID[maxdim];
  Int_t MCTrackID[maxdim];
  Int_t MCMotherID[maxdim];
/*
  Int_t   trid=0;
  Int_t   nseg=0;
  Int_t   npl=0;
  Int_t   n0=0;
  Float_t xv=0.;
  Float_t yv=0.;*/

  int nentr = (int)(vtx->GetEntries());
  //if(cut) Log(2,"EdbDataProc::ReadVtxTree","select %d of %d vertices by cut %s",nlst, nentr, cut.GetTitle() );

  TClonesArray *tracks = new TClonesArray("EdbSegP");
  TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
  TClonesArray *segf = new TClonesArray("EdbSegP", 60);
  
  EdbSegP *trk=0;
  EdbSegP *s1=0;
  //vertex variables
  vtx->SetBranchAddress("vID",&vID);
  vtx->SetBranchAddress("flag",&flag);
  vtx->SetBranchAddress("vx",&vx);
  vtx->SetBranchAddress("vy",&vy);
  vtx->SetBranchAddress("vz",&vz);
  vtx->SetBranchAddress("n",&n);
  //arrays
  vtx->SetBranchAddress("nseg",&nseg);
  vtx->SetBranchAddress("TrackID",&TrackID);
  vtx->SetBranchAddress("incoming",&incoming);
  //TClones arrays
  vtx->SetBranchAddress("sf", &segf);
  vtx->SetBranchAddress("s",  &seg);
  vtx->SetBranchAddress("t.", &tracks);

  EdbPattern *pat=0;
  int entr=0;
  //now each tree entry is a vertex

  TObjArray * vertextracks = new TObjArray(maxdim);
  float expz = 0.;
  EdbVertexRec *vertexrec = new EdbVertexRec(); //test to use addtrack method
  vertexrec->SetPVRec(&ali);
  //RESET COUNTERS!
  int itotalseg = 0; //counter for all the segments
  vertextracks->Clear("C");
  vtx->BuildIndex("vID");
  vtx->GetEntryWithIndex(vertexID);
  EdbVertex *v1 = new EdbVertex();
  //loop on all tracks associated to the vertex
  for (int itrk = 0; itrk < n; itrk++){
     EdbTrackP *tr1 = new EdbTrackP();
     trk = (EdbSegP*) tracks->At(itrk); //getting track itrk
     ((EdbSegP*)tr1)->Copy(*trk);
     tr1->SetM(0.139);                 //TODO

     for(int i=0; i<nseg[itrk]; i++) {
      s1 = (EdbSegP*)(seg->At(itotalseg));
      pat = ali.GetPattern( s1->PID() );
      if(!pat) { 
	     //Log(1,"EdbDataProc::ReadTracksTree","WARNING: no pattern with pid %d: creating new one!",s1->PID()); 
	     pat = new EdbPattern( 0., 0., s1->Z() );
	     pat->SetID(s1->PID());
	     pat->SetScanID(s1->ScanID());
	     ali.AddPatternAt(pat,s1->PID());
      }

       tr1->AddSegment( pat->AddSegment(*s1 ) );
       tr1->AddSegmentF( new EdbSegP(*((EdbSegP*)(segf->At(itotalseg)))) );
       itotalseg++; //ALWAYS REMEMBER COUNTERS, YOU IDIOT!
    }
     tr1->SetSegmentsTrack(tr1->ID());
     tr1->SetCounters();
    //tr1->FitTrackKFS(true);
     tr1->SetTrack(TrackID[itrk]); //providing trackid to eTrack so it will not be lost when trackID resets
     ali.AddTrack(tr1);
     vertextracks->Add(tr1);     
     //v1 = vertexrec->AddTrackToVertex(v1, tr1, incoming[itrk]);
    }
    //setting vertex parameters and saving vertex
  v1 = vertexrec->Make1Vertex(*vertextracks,vz);
  v1->SetID(vID);
  v1->SetFlag(flag);
  v1->SetXYZ(vx,vy,vz);   
  ali.AddVertex(v1);
  return v1;
//  Log(2,"EdbDataProc::ReadVtxTree","%d tracks are read",nlst);
 }
