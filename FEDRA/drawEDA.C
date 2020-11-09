EdbVertex* GetVertexFromTree( EdbPVRec &ali, const char     *fname, const int vertexID );

EdbPVRec * dummyrec = new EdbPVRec();
EdbEDA * eda = new EdbEDA(dummyrec,false); // init DataSet but doesn't read linked_track.root
//EdbEDA *eda;

void loadcouples(EdbPVRec * ali, float xcenter, float ycenter, float rmax = 2000);
void showerEDA(const int trackID, bool first = true, TString trackfilename = "linked_tracks.root");
void showerall(){
 showerEDA(1388);
 eda->Reset();
 showerEDA(99188);
 //showerEDA(99188,false);
}

void showerEDA(const int trackID, bool first = true, TString trackfilename = "linked_tracks.root"){

 bool drawtracks = true;
 const int ntracks = 1;
 int tracklist[ntracks] = {trackID}; //87251
 
 EdbPVRec *ali = new EdbPVRec();

 EdbDataProc* dproc = new EdbDataProc();

 dproc->InitVolume(100, "nseg>1 && trid<20000"); //reading some tracks to set the patterns
 ali = dproc->PVR();
 ali->eTracks=0; //removing tracks used for patterns
 ali->FillCell(30,30,0.009,0.009);
float xcenter, ycenter; //shower center
int startplate;
EdbSegP *selectorsegment; //source of shower
  //reading track file and setting branches
if (drawtracks){
 TFile *trackfile = TFile::Open(trackfilename.Data());
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
 EdbSegP *s1f=0;

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
 EdbPattern *pat = 0;
 for (int trackID: tracklist){
  //if (!(ali->GetTrack(trackID))) continue;
  seg->Clear();
  segf->Clear();
  tracks->GetEntryWithIndex(trackID);

  EdbTrackP *tr1 = new EdbTrackP();
  ((EdbSegP*)tr1)->Copy(*trk);
  tr1->SetM(0.139);                 //TODO  
    for(int i=0; i<nseg; i++) {
      s1 = (EdbSegP*)(seg->At(i));
      pat = ali->GetPattern(s1->PID());
      if (!pat){
        cout<<"Creating new pattern"<<endl;
        pat = new EdbPattern( 0., 0., s1->Z());
        pat->SetID(s1->PID());
        pat->SetScanID(s1->ScanID());
        ali->AddPatternAt(pat,s1->PID());
      }

      s1f = (EdbSegP*)(segf->At(i));
      s1->SetDZ(300);
      s1f->SetDZ(300);
      s1->SetFlag(kMagenta); //EDA colours tracks according to flag
      tr1->AddSegment(  new EdbSegP(*((EdbSegP*)(s1))) );
      tr1->AddSegmentF( new EdbSegP(*((EdbSegP*)(s1f))) );
    } //end on loop on segments associated to the track
  tr1->SetSegmentsTrack(tr1->ID());
  tr1->SetCounters();    

  //adding track to drawing set and building EDA
  ali->AddTrack(tr1);

  //setting couples region
  xcenter = tr1->GetSegmentFirst()->X();
  ycenter = tr1->GetSegmentFirst()->Y();
  startplate = tr1->GetSegmentFirst()->Plate();
  } //end loop on tracks
 }
 //test loading couples
 loadcouples(ali, xcenter,ycenter,5000);
 eda->SetPVR(ali);
 eda->GetTrackSet("TS")->AddTracksPVR(ali);
 //eda = new EdbEDA(ali);
 EdbEDAUtil::FillTracksFromPatterns(ali);
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
 eda->Run();

 TObjArray * segarray = new TObjArray(100);
 segarray->Add(eda->GetTrackSet("TS")->GetTrack(0)->GetSegmentFirst());
 eda->SetSelected(segarray);

 //test showering by command
 EdbEDAShowerTab *showtab = new EdbEDAShowerTab();
 showtab->Reco();

 TFile* file=TFile::Open("Shower.root");
 TTree* ShowerTree = (TTree *)file->Get("treebranch");

 int showersize;
 float output15, output30;

 ShowerTree->SetBranchAddress("sizeb",&showersize);
  ShowerTree->SetBranchAddress("output15",&output15);
 ShowerTree->SetBranchAddress("output30",&output30);
 int nshowers = ShowerTree->GetEntries();
 cout<<"Shower Selector: Segment from plate"<<startplate<<endl;
 if (nshowers == 0) cout<<"No shower found!"<<endl;
 else{
   cout<<"Startplate\tShower size\toutput15\toutput30"<<endl;
   for (int ishower = 0; ishower < nshowers; ishower++){
     ShowerTree->GetEntry(ishower);
     cout<<startplate<<"\t"<<showersize<<"\t"<<output15<<"\t"<<output30<<endl;
   }



 }
 //delete eda;
 //delete ali;
}

void drawEDA(bool newversion = true, TString vertexfilename = " vertextree.root", TString trackfilename = "linked_tracks.root"){
 bool drawvertices = true;
 bool drawtracks = false;
 const int nvertices = 1;
 int vertexlist[nvertices] = {110827};
 int vertexcolors[nvertices] = {kBlue};
 const int ntracks = 1;
 int tracklist[ntracks] = {99188}; //87251
 
 EdbPVRec *ali = new EdbPVRec();

 EdbDataProc* dproc = new EdbDataProc();

 dproc->InitVolume(100, "nseg>1 && trid<20000"); //reading some tracks to set the patterns
 ali = dproc->PVR();
 ali->eTracks=0; //removing tracks used for patterns
 ali->FillCell(30,30,0.009,0.009);
 cout<<"TEST"<<endl;
 EdbVertex *vertex = NULL;
 if (drawvertices){
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int vID = vertexlist[ivtx];
 //f or (int vID: vertexlist){ //range for loop, C++11
    if (newversion) vertex = dproc->GetVertexFromTree(*ali,vertexfilename,vID);
    else{ 
      TFile * inputfile = TFile::Open(vertexfilename.Data());
      EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");
      vertex = (EdbVertex*) vertexrec->eVTX->At(vID);
    }
    for (int itrk = 0; itrk < vertex->N(); itrk++){
      EdbTrackP* track =  vertex->GetTrack(itrk);
      bool foundtrack = false;
      for (int trackID: tracklist){      
 
       if (track->Track()==trackID){
          foundtrack = true;
          for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(kMagenta); // to color them differently    }
       } 
      if (!foundtrack)
        for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[ivtx]); // to color them differently
      }
    }
  }
 }

float xcenter, ycenter; //shower cente
EdbSegP *selectorsegment; //source of shower
  //reading track file and setting branches
if (drawtracks){
 TFile *trackfile = TFile::Open(trackfilename.Data());
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
 EdbSegP *s1f=0;

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
 EdbPattern *pat = 0;
 for (int trackID: tracklist){
  //if (!(ali->GetTrack(trackID))) continue;
  seg->Clear();
  segf->Clear();
  tracks->GetEntryWithIndex(trackID);

  EdbTrackP *tr1 = new EdbTrackP();
  ((EdbSegP*)tr1)->Copy(*trk);
  tr1->SetM(0.139);                 //TODO  
    for(int i=0; i<nseg; i++) {
      s1 = (EdbSegP*)(seg->At(i));
      pat = ali->GetPattern(s1->PID());
      if (!pat){
        cout<<"Creating new pattern"<<endl;
        pat = new EdbPattern( 0., 0., s1->Z());
        pat->SetID(s1->PID());
        pat->SetScanID(s1->ScanID());
        ali->AddPatternAt(pat,s1->PID());
      }

      s1f = (EdbSegP*)(segf->At(i));
      s1->SetDZ(300);
      s1f->SetDZ(300);
      s1->SetFlag(kMagenta); //EDA colours tracks according to flag
      tr1->AddSegment(  new EdbSegP(*((EdbSegP*)(s1))) );
      tr1->AddSegmentF( new EdbSegP(*((EdbSegP*)(s1f))) );
    } //end on loop on segments associated to the track
  tr1->SetSegmentsTrack(tr1->ID());
  tr1->SetCounters();    

  //adding track to drawing set and building EDA
  ali->AddTrack(tr1);

  //setting couples region
  xcenter = tr1->GetSegmentFirst()->X();
  ycenter = tr1->GetSegmentFirst()->Y();
  } //end loop on tracks
 }
 //test loading couples
 loadcouples(ali, xcenter,ycenter);
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 EdbEDAUtil::FillTracksFromPatterns(ali);
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
 eda->Run();

 TObjArray * segarray = new TObjArray(100);
 segarray->Add(eda->GetTrackSet("TS")->GetTrack(0)->GetSegmentFirst());
 eda->SetSelected(segarray);
/*
 //test showering by command
 EdbEDAShowerTab *showtab = new EdbEDAShowerTab();
 showtab->Reco();

 TFile* file=TFile::Open("Shower.root");
 TTree* ShowerTree = (TTree *)file->Get("treebranch");

 int showersize;
 float output15, output30;

 ShowerTree->SetBranchAddress("sizeb",&showersize);
  ShowerTree->SetBranchAddress("output15",&output15);
 ShowerTree->SetBranchAddress("output30",&output30);
 int nshowers = ShowerTree->GetEntries();
 if (nshowers == 0) cout<<"No shower found!"<<endl;
 else{
   cout<<"Shower size \t output30"<<endl;
   for (int ishower = 0; ishower < nshowers; ishower++){
     ShowerTree->GetEntry(ishower);
     cout<<showersize<<"\t"<<output15<<"\t"<<output30<<endl;
   }


 }*/
/* TFile *file = new TFile("File_EdbPVRec.root","RECREATE");
 ali->Write();
 file->Close()*/;

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
       itotalseg++;
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


void loadcouples(EdbPVRec * ali, float xcenter, float ycenter, float rmax = 2000){
  //loading couples at a distance of rmax between xcenter and ycenter
  //float xcenter = 76145.;
  //float ycenter = 33233.;
  TString runpath = TString("./");

  //float rmax = 2000; //maximum distance in xy

  TFile *setfile = TFile::Open("b000001.0.0.0.set.root");
  EdbScanSet *set = (EdbScanSet*) setfile->Get("set");
//  TString condition = TString("1");
 TString condition = TString(Form("TMath::Sqrt(pow(s.eX - %f,2)+ pow(s.eY - %f,2))<%f",xcenter,ycenter,rmax));
 //loop over plates
  const int nplates=29;
  EdbCouplesTree *ect[nplates];
  for (int i = nplates; i >= 1; i--){ 
  //getting z position and affine transformation
  EdbPlateP* p = set->GetPlate(i);
  float zplate = p->Z();
 /* //if(!p) continue;
  EdbAffine2D *aff = new EdbAffine2D();
  set->GetAffP2P(i, nplates, *aff); //usually last plate is the reference one
  aff->Print();*/ //affine transformations
  ect[i-1] = new EdbCouplesTree();
  //getting couples
  if (i <10) ect[i-1]->InitCouplesTree("couples",(runpath+TString(Form("couples/p00%i/1.%i.0.0.cp.root",i,i))).Data(),"READ");
  else ect[i-1]->InitCouplesTree("couples",(runpath+TString(Form("couples/p0%i/1.%i.0.0.cp.root",i,i))).Data(),"READ");

  //loop into couples (only the ones passing condition)
  ect[i-1]->eTree->Draw(">>goodcouples", condition.Data());
  TEventList *goodcouples = (TEventList*) gDirectory->GetList()->FindObject("goodcouples");
  
  const int ngoodcouples = ect[i-1]->eTree->GetEntries(condition.Data());

  //int nseg = ect[i-1]->eTree->GetEntries();
  cout<<"Reading "<<ngoodcouples<<" from plate "<<i<<endl;
  //*************LOOP OVER SEGMENTS***************
  for (int iseg = 0; iseg< ngoodcouples; iseg++){
   //getting entry of good segment;
   int igoodsegment = goodcouples->GetEntry(iseg);
   //***Getting information about that segment***;
   ect[i-1]->GetEntry(igoodsegment);
   EdbSegP *seg = new EdbSegP();
   seg->Copy(*(ect[i-1]->eS));
   seg->SetZ(zplate);
   seg->SetPlate( i ); //setting which plate this couple belongs to
   seg->SetPID(29-i);
   ali->AddSegment(*seg);
    } //end loop over segments
  }//end loop over plates
 for (int ipattern=1; ipattern<ali->Npatterns(); ipattern++) {
   ali->GetPattern(ipattern)->SetPID(29-ipattern); 
        }
}