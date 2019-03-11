//----------------------------------------------------------------------------
//
//  Usage: $ root -l 
//         root[1] .L check_vertex.C
//         root[2] trseg(16)  // process the event 16 (simulation of cause!)
//
//  Usage: $ root -l 
//         root[1] .L check_vertex.C
//         root[2] trvol()  // do all vertex reconstruction starting from the linked_tracks.root
//
//----------------------------------------------------------------------------

EdbDataProc  *dproc=0;
EdbPVRec     *gAli=0;
EdbVertexRec *gEVR=0;
EdbDisplay   *ds=0;

TH1D *hip = new TH1D("hip","Impact parameters",500,0,5000);

void trvol( const char *def, const char *rcut = "nseg>1" );
void init( const char *def, int iopt,  const char *rcut="1" );
void set_segments_dz(float dz=300.);
void do_propagation();
void do_vertex();
void td();
void sd();
void vd( int trmin, float amin);

#include <vector>
namespace VERTEX_PAR
{
  float DZmax = 3000.;
  //float ProbMinV   = 0.0001;  // minimum acceptable probability for chi2-distance between tracks
  float ProbMinV   = 0.01;
  float ImpMax     = 15.;    // maximal acceptable impact parameter [microns] (for preliminary check)
  bool  UseMom     = false;  // use or not track momentum for vertex calculations
  bool  UseSegPar  = true;  // use only the nearest measured segments for vertex fit (as Neuchatel)
  int   QualityMode= 0;      // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
}

//---------------------------------------------------------------------
void charm_vertexing(char *dset=0)
{
  //trvol(dset);
  //trvol(dset,"nseg>1 && TMath::Abs(s.eY-50000)<5000 && TMath::Abs(s.eX-50000)<5000");   
// SELECTION WHICH TRACKS TO USE FOR VERTEXING
 trvol(dset,"nseg>1 && TMath::Abs(s.eY-50000)<5000");

   // reconstruct vertexes starting from linked_tracks.root
}

//---------------------------------------------------------------------
void trvol( const char *def, const char *rcut )
{
  // this function read volume tracks and do the vertex reconstruction
  // from linked_tracks.root

  init(def, 100 ,rcut);                      // read tracks (option 100)
  gAli->FillCell(30,30,0.009,0.009);
  //do_propagation();
  do_vertex();
  //vd(4,0.01);   // draw reconstructed vertex 
  //td(); //to draw tracks
  //dproc->MakeTracksTree(gAli,"linked_tracks_p.root");
}

//---------------------------------------------------------------------
void init( const char *def, int iopt,  const char *rcut)
{
  if(!def)  dproc = new EdbDataProc();
  else      dproc = new EdbDataProc(def);

  dproc->InitVolume(iopt, rcut);  
  gAli = dproc->PVR();
  set_segments_dz(300.);
}

//---------------------------------------------------------------------
void set_segments_dz(float dz)
{
  int np = gAli->Npatterns();
  for(int i=0; i<np; i++) {
    EdbPattern *p = gAli->GetPattern(i);
    int ns = p->N();
    for(int j=0; j<ns; j++) p->GetSegment(j)->SetDZ(dz);
  }
}

//---------------------------------------------------------------------
void do_propagation()
{
  // example of additional propagation and other tracking operations if necessary
  EdbTrackFitter tf;

  EdbTrackP *tr=0;
  int ntr = gAli->eTracks->GetEntries();
  printf("ntr = %d\n",ntr);

  for(int i=0; i<ntr; i++) {
    tr = (EdbTrackP*)(gAli->eTracks->At(i));
    tr->SetID(i);
    tr->SetSegmentsTrack();
  }

  int nadd = 0;
  int nseg=0;
  
  for(int i=0; i<ntr; i++) {
    tr = (EdbTrackP*)(gAli->eTracks->At(i));
    
    //float p=momentum; //not working on ROOT6
    //float p=tf.P_MS(*tr);
    //if(tr->ID()==176) p=0.33;
    //if(tr->ID()==230) p=22.;    

    //tr = tr->SetErrorP(0.2*0.2*p*p); //not working on ROOT6
    nseg = tr->N();
    //tr->SetP(p); //for the previous one
    if(tr->Flag()<0) continue;
    nadd += gAli->PropagateTrack( *tr, true, 0.01, 3, 0 );
    if(tr->Flag()<0) printf("%d flag=%d\n",i,tr->Flag());
    //if(tr->N() != nseg) printf("%d nseg=%d (%d) \t p = %f\n",i,tr->N(),nseg,tr->P());
  }
  printf("nadd = %d\n",nadd);
}

//---------------------------------------------------------------------
void do_vertex()
{
  using namespace VERTEX_PAR;

  //gAli->FitTracks(4.,0.139 );
  gEVR = new EdbVertexRec();
  gEVR->eEdbTracks = gAli->eTracks;
  gEVR->eVTX       = gAli->eVTX;
  gEVR->SetPVRec(gAli);

  gEVR->eDZmax=DZmax;
  gEVR->eProbMin=ProbMinV;
  gEVR->eImpMax=ImpMax;
  gEVR->eUseMom=UseMom;
  gEVR->eUseSegPar=UseSegPar;
  gEVR->eQualityMode=QualityMode;

  printf("%d tracks vor vertexing\n",  gEVR->eEdbTracks->GetEntries() );
  int nvtx = gEVR->FindVertex();
  printf("%d 2-track vertexes was found\n",nvtx);

  if(nvtx == 0) return;
  int nadd =  gEVR->ProbVertexN();

  int nl = gEVR->LinkedVertexes(); //should avoid same track associated to multiple vertices
  printf("%d vertices are linked\n",nl);

  TFile *fvtx = new TFile("vertices_MC.root","RECREATE"); //OUTPUT FILE NAME
  TTree *vtx = new TTree("vtx","vtx");
  Float_t vx, vy, vz;
  Int_t vID = 0;
  Float_t maxaperture;
  Float_t probability;
  Int_t n;
  const Int_t maxdim = 1000;
  //big arrays for containers of track variables
  Int_t nholes[maxdim];
  Int_t maxgap[maxdim];
  Int_t nseg[maxdim];
  Float_t TX[maxdim];
  Float_t TY[maxdim];
  Float_t impactparameter[maxdim];
  Int_t incoming[maxdim]; //coming from the vertex
 
  vtx->Branch("vID",&vID,"vID/I");
  vtx->Branch("vx",&vx,"vx/F");
  vtx->Branch("vy",&vy,"vy/F");
  vtx->Branch("vz",&vz,"vz/F");
  vtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
  //vtx->Branch("maxrmsthetaspace",&maxrmsthetaspace,"maxrmsthetaspace/F");
  vtx->Branch("probability",&probability,"probability/F");
  vtx->Branch("n",&n,"n/I");
  //track variables (they are array with the number of tracks as size)
  vtx->Branch("TX",&TX,"TX[n]/F");
  vtx->Branch("TY",&TY,"TY[n]/F");
  vtx->Branch("nseg",&nseg,"nseg[n]/I");
  vtx->Branch("nholes",&nholes,"nholes[n]/I"); //even more variables in this tree
  vtx->Branch("maxgap",&maxgap,"maxgap[n]/I"); 
  vtx->Branch("incoming",&incoming,"incoming[n]/I");
  vtx->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
  
  EdbVertex *vertex = new EdbVertex();
  EdbTrackP *track = new EdbTrackP();

  //gAli->Write();
  fvtx->cd();
  TObjArray *varr = new TObjArray(); //vertices to be saved
  //gEVR->Write();
  cout<<"Ho salvato i vertici nel file"<<endl;
  for(Int_t ivtx=0; ivtx<gEVR->eVTX->GetEntries(); ivtx++){
    vertex=(EdbVertex*)(gEVR->eVTX->At(ivtx));
    vx=vertex->X();
    vy=vertex->Y();
    vz=vertex->Z();
    n=vertex->N();
    maxaperture = vertex->MaxAperture();
    probability = vertex->V()->prob();  
    if(vertex->Flag()<0) continue; //saving only 'true' vertices in the tree file and in the object
    if(n < 4) continue; 
    if (maxaperture <0.01) continue;
    varr->Add(vertex); //true vertex, we can save it
    //adding vertex to list to be saved
    //loop on tracks //now it can be done offline (again)
    /*for (int itrk = 0; itrk < n; itrk++){
     track = vertex->GetTrack(itrk);
     nseg[itrk] = track->N();
     Int_t zpos = vertex->GetVTa(itrk)->Zpos();
     incoming[itrk] = zpos;
     TX[itrk] = track->TX();
     TY[itrk] = track->TY();
     nholes[itrk] = track->N0();
     maxgap[itrk] = track->CheckMaxGap();
     impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
    }*/
    vtx->Fill();
    vID++; //now the number of elements in the tree and vertices is the same. vID starts from 0
  }
  //trying to save the selected vertices
  EdbVertexRec *mygEVR = new EdbVertexRec();
  mygEVR->eVTX = varr;  
  //using same parameters for vertexrec as the original one
  mygEVR->eDZmax=DZmax;
  mygEVR->eProbMin=ProbMinV;
  mygEVR->eImpMax=ImpMax;
  mygEVR->eUseMom=UseMom;
  mygEVR->eUseSegPar=UseSegPar;
  mygEVR->eQualityMode=QualityMode;
  mygEVR->Write();

  vtx->Write();
  hip->Draw();
  fvtx->Close(); //close the file where vertices are saved
  //fclose(tvtx);
}

//---------------------------------------------------------------------
void td()
{
  // draw all tracks

  TObjArray *trarr=gAli->eTracks;
  gStyle->SetPalette(1);
  const char *dsname="display-t";
  ds = EdbDisplay::EdbDisplayExist(dsname);
  if(!ds)  ds=new EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.);
  ds->SetVerRec(gEVR);
  ds->SetDrawTracks(4);
  ds->SetArrTr( trarr );
  printf("%d tracks to display\n", trarr->GetEntries() );
  ds->Draw();
}

//---------------------------------------------------------------------
void sd()
{
  // draw all tracks and segments (basetracks)

  TObjArray *trarr = gAli->eTracks;
  TObjArray *sarr  = new TObjArray();

  EdbSegP *s=0;
  for(int i=0; i<gAli->Npatterns(); i++) {
    EdbPattern *pat = gAli->GetPattern(i);
    for(int j=0; j<pat->N(); j++) {
      s = pat->GetSegment(j);
      if(s->Track()<0)                  // exclude segments already attached to tracks
	sarr->Add(s);
    }
  }

  printf("%d tracks to display\n",   trarr->GetEntries() );
  printf("%d segments to display\n", sarr->GetEntries()  );

  gStyle->SetPalette(1);
  const char *dsname="display-s";
  ds = EdbDisplay::EdbDisplayExist(dsname);
  if(!ds)  ds=new EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.);
  ds->SetVerRec(gEVR);
  ds->SetDrawTracks(4);
  ds->SetArrSegP( sarr );
  ds->SetArrTr( trarr );
  ds->Draw();
}

//---------------------------------------------------------------------
void vd( int trmin, float amin)
{
  // draw vertexes with multiplicity>=trmin, and aperture >= amin

  TObjArray *varr = new TObjArray();
  TObjArray *tarr = new TObjArray();

  EdbVertex *v=0;
  EdbTrackP *t=0;

  int nv = gEVR->Nvtx();
  printf("nv=%d\n",nv);
  if(nv<1) return;

  for(int i=0; i<nv; i++) {
    v = (EdbVertex *)(gEVR->eVTX->At(i));
    if(v->Flag()<0)         continue;
    if( v->N()<trmin )                       continue;
    if( v->N()<3 && v->MaxAperture()<amin )  continue;

    varr->Add(v);
    for(int j=0; j<v->N(); j++) tarr->Add( v->GetTrack(j) );
  }

  gStyle->SetPalette(1);

  const char *dsname="display-v";
  ds = EdbDisplay::EdbDisplayExist(dsname);
  if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-40000., 0.);
  ds->SetVerRec(gEVR);
  ds->SetArrTr( tarr );
  printf("%d tracks to display\n", tarr->GetEntries() );
  ds->SetArrV( varr );
  printf("%d vertex to display\n", varr->GetEntries() );
  ds->SetDrawTracks(4);
  ds->SetDrawVertex(1);
  ds->Draw();
}

// Experimental function to study segment distance to vertex
void add_proton(EdbVertex* vertex, EdbVertexRec* gEVR){
cout<<"start analyzing a vertex"<<endl;
   float track_ipmin = 10000.;
   int ntr  = gEVR->eEdbTracks->GetEntries(); //number of tracks 
//   EdbVTA *vta = 0;
   EdbSegP *potsegment;
   
   bool associated_track; //track near the vertex

   for(int itr=0; itr<ntr; itr++)   {
     associated_track = false;
     EdbTrackP *tr = (EdbTrackP*)(gEVR->eEdbTracks->At(itr));
    
     if (tr->N() > 26){

	float iptrack2vertex = vertex->DistTrack(tr,0.);
	if (iptrack2vertex < track_ipmin) track_ipmin = iptrack2vertex;

	if (iptrack2vertex < 200.) associated_track = true; //this is not the right track

	if (associated_track) cout<<"Let's see each segment for this track, track distance = "<<iptrack2vertex<<endl;

	for (int iseg = 0; iseg < tr->N(); iseg++){ //loop on segments

	   potsegment = (EdbSegP*)tr->GetSegment(iseg);
	   float ipseg2vertex = vertex->DistSeg(potsegment,0.);         
	   if (associated_track) cout<<ipseg2vertex<<" "<<potsegment->Z()<<endl;         	           
	}	
     }   
 if (associated_track) cout<<endl;
   } 
 hip->Fill(track_ipmin);
}
