#include "TObject.h"
#include "EdbPVRec.h"
#include "EdbPattern.h"
#include "EdbVertex.h"
#include "TFile.h"
#include "TCut.h"
#include "TEventList.h"

// structure of vertextree
/* 
 Each entry is a vertex
 float and int numbers for vertex variables
 TClonesArrays of EdbSegP objects for tracks, segments and fitted segments
 
 Size of t. array:
  number of tracks in the vertex

 Size of s and sf arrays:
  total number of segments and fitted segments (each track has many segments) 

*/
class VertexIO:public TObject{
 public:
 static int MakeVertexTree(TObjArray &vtxarr, const char *file)
 {
  //Log(2,"MakeVertexTree","write vertices into %s ... ",file);
  TFile fil(file,"RECREATE");
  TTree *vtx= new TTree("vtx","Reconstructed vertices in emulion");
  
  EdbSegP      *tr = new EdbSegP();
  EdbTrackP *track = NULL;

  TClonesArray *tracks = new TClonesArray("EdbSegP"); //tracks are saved as EdbSegP, like in MakeTracksTree
  TClonesArray *segments  = new TClonesArray("EdbSegP");
  TClonesArray *segmentsf = new TClonesArray("EdbSegP");

  //list of variables

  Float_t vx, vy, vz;
  Int_t vID = 0;
  Float_t maxaperture;
  Float_t probability;
  Int_t n;
  Int_t flag; //we need also the flag to check vertex type
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

  //list of branches
  vtx->Branch("vID",&vID,"vID/I");
  vtx->Branch("flag",&flag,"flag/I");
  vtx->Branch("vx",&vx,"vx/F");
  vtx->Branch("vy",&vy,"vy/F");
  vtx->Branch("vz",&vz,"vz/F");
  vtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
  //vtx->Branch("maxrmsthetaspace",&maxrmsthetaspace,"maxrmsthetaspace/F");
  vtx->Branch("probability",&probability,"probability/F");
  vtx->Branch("n",&n,"n/I");
  vtx->Branch("t.",&tracks); // track data is now filled after using Copy method
  //tracks->Branch("t.","EdbSegP",&track,32000,99);
  vtx->Branch("s", &segments);
  vtx->Branch("sf",&segmentsf);
  //track variables (they are array with the number of tracks as size)
  vtx->Branch("TrackID",&TrackID,"TrackID[n]/I");
  vtx->Branch("nseg",&nseg,"nseg[n]/I");
  vtx->Branch("nholes",&nholes,"nholes[n]/I"); //even more variables in this tree
  vtx->Branch("maxgap",&maxgap,"maxgap[n]/I"); 
  vtx->Branch("incoming",&incoming,"incoming[n]/I");
  vtx->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
  //inserting MCtrue information
  vtx->Branch("MCEventID", &MCEventID, "MCEventID[n]/I");
  vtx->Branch("MCTrackID",&MCTrackID,"MCTrackID[n]/I");
  vtx->Branch("MCMotherID",&MCMotherID,"MCMotherID[n]/I");

  int nvtx = vtxarr.GetEntriesFast();
  //START TREE FILLING
  for(int ivtx=0; ivtx<nvtx; ivtx++) {
    int itotalseg = 0;
    tracks->Clear("C");
    segments->Clear("C");
    segmentsf->Clear("C");

    EdbVertex* vertex = (EdbVertex*) (vtxarr.At(ivtx));
    if(vertex->Flag()<0) continue; //saving only 'true' vertices in the tree file and in the object
    vx=vertex->X();
    vy=vertex->Y();
    vz=vertex->Z();
    n=vertex->N();
    maxaperture = vertex->MaxAperture();
    probability = vertex->V()->prob();  
    flag = vertex->Flag();
    //adding vertex to list to be saved
    //loop on tracks //now it can be done offline (again)
    for (int itrk = 0; itrk < n; itrk++){
     //getting track and track variables to fill branches
     track = vertex->GetTrack(itrk);
     tr->Copy(*track);
     if(tr) new((*tracks)[itrk])  EdbSegP( *tr ); //adding track to trackclonesarray

     TrackID[itrk] = track->Track(); //the eTrack attribute of EdbSegP now allows association to original root tree
     nseg[itrk] = track->N();
     Int_t zpos = vertex->GetVTa(itrk)->Zpos();
     incoming[itrk] = zpos;
     nholes[itrk] = track->N0();
     maxgap[itrk] = track->CheckMaxGap();
     impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
     //Storing MCTrue information (of course for real data these values have no sense)
     MCEventID[itrk] = track->MCEvt();
     MCTrackID[itrk] = track->MCTrack();
     MCMotherID[itrk] = track->Aid(0); //used to store MotherID information
    //w    = track->Wgrains();
     EdbSegP *s=0,*sf=0;
    //loop on segments
     for(int is=0; is<nseg[itrk]; is++) {
      s = track->GetSegment(is);     
      if(s) new((*segments)[itotalseg])  EdbSegP( *s );
      sf = track->GetSegmentF(is);
      if(sf) new((*segmentsf)[itotalseg])  EdbSegP( *sf );
      itotalseg++;
      }

   // track->SetVid( 0, tracks->GetEntries() );  // put track counter in t.eVid[1]
     }
    vtx->Fill();
    vID++; //now the number of elements in the tree and vertices is the same. vID starts from 0
    }
  vtx->Write();
  fil.Close();
  //Log(2,"EdbDataProc::MakeVertexTree","%d vertices are written",nvtx);
  return nvtx; 
}


static int ReadVertexTree( EdbPVRec &ali, const char     *fname, const char *rcut )
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
  TCut cut = rcut;
  vtx->Draw(">>lst", cut );
  TEventList *lst = (TEventList*)gDirectory->GetList()->FindObject("lst");
  int nlst =lst->GetN();
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
  for (int j=0; j<nlst; j++){
    //RESET COUNTERS!
    int itotalseg = 0; //counter for all the segments
    vertextracks->Clear("C");

    entr = lst->GetEntry(j);
    vtx->GetEntry(entr);
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
  }

//  Log(2,"EdbDataProc::ReadVtxTree","%d tracks are read",nlst);
  return nlst;
 }

 static EdbVertex* GetVertexFromTree( EdbPVRec &ali, const char     *fname, const int vertexID )
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
 ClassDef(VertexIO,1);
};
 ClassImp(VertexIO)
