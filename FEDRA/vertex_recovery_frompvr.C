//trying to recover vT vertex: EdbVertexRec::MakeV( EdbVertex &edbv, bool isRefit )

//values for the fit as when the vertexing was performed

namespace VERTEX_PAR
{
  float DZmax      = 3000.;  // maximum z-gap in the track-vertex group
  float ProbMinV   = 0.01;  // minimum acceptable probability for chi2-distance between tracks
  float ImpMax     = 15.;    // maximal acceptable impact parameter [microns] (for preliminary check)
  bool  UseMom     = false;  // use or not track momentum for vertex calculations
  bool  UseSegPar  = false;  // use only the nearest measured segments for vertex fit (as Neuchatel)
  int   QualityMode= 0;      // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
}

EdbDataProc *dproc = new EdbDataProc();
EdbPVRec *gAli = dproc->PVR();
void vertex_recovery_frompvr(){
 TFile *vertexfile = TFile::Open("/ship/CHARM2018/CH6-R2/b000001/vertices_uppertracks.root");
 EdbPVRec *pvrec = vertexfile->Get("EdbPVRec");
 TObjArray *vertexlist = pvrec->eVTX;
 //values for the fit as when the vertexing was performed
 using namespace VERTEX_PAR;
 EdbVertexRec *vertexrec = new EdbVertexRec();
 gEVR->eDZmax=DZmax;
 gEVR->eProbMin=ProbMinV;
 gEVR->eImpMax=ImpMax;
 gEVR->eUseMom=UseMom;
 gEVR->eUseSegPar=UseSegPar;
 gEVR->eQualityMode=QualityMode;

 vertexrec->eVTX       = vertexlist;

 vertexrec->SetPVRec(pvrec);
 //arrays for vertices drawing
 TObjArray *varr = new TObjArray();
 TObjArray *tarr = new TObjArray();
 int nvertices = vertexrec->Nvtx();

 //main loop
 for (int i = 0; i < nvertices; i++){
  EdbVertex *vertex = vertexrec->GetVertex(i);
  vertexrec->MakeV(*vertex);
  varr->Add(vertex);
  for (int j= 0; j<vertex->N();j++){
   EdbTrackP * track = vertex->GetTrack(j);
   tarr->Add(track);
   
  }
 }

 //drawing the results
 const char *dsname="display-v";
 EdbDisplay   *ds=0;
 ds = EdbDisplay::EdbDisplayExist(dsname);
 if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-40000., 0.);
 ds->SetVerRec(vertexrec);
 ds->SetArrTr( tarr );
 printf("%d tracks to display\n", tarr->GetEntries() );
 ds->SetArrV( varr );
 printf("%d vertex to display\n", varr->GetEntries() );
 ds->SetDrawTracks(4);
 ds->SetDrawVertex(1);
 ds->Draw();
}
