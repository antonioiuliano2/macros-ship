//display for vertices, patterns are not drawn
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

void drawEDAvertex(bool newversion = true, TString vertexfilename= "vertextree_test.root"){
 using namespace VERTEX_PAR;
 const int nvertices = 1;
 int vertexlist[nvertices] = {7197};
 int vertexcolors[nvertices] = {kRed};
 EdbDataProc *dproc = new EdbDataProc();

 TFile * inputfile = TFile::Open(vertexfilename.Data());
 EdbVertexRec *vertexrec; 

 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);

// EdbTrackP *specialtrack = new EdbTrackP();

 EdbPVRec *ali = new EdbPVRec();
 EdbScanCond *scancond = new EdbScanCond();
 ali->SetScanCond(scancond);

 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = vertexlist[i];
  
  EdbVertex *vertex = 0;

  if (newversion){ 
    vertexrec = new EdbVertexRec();
    vertexrec->SetPVRec(ali);
    vertexrec->eDZmax=DZmax;
    vertexrec->eProbMin=ProbMinV;
    vertexrec->eImpMax=ImpMax;
    vertexrec->eUseMom=UseMom;
    vertexrec->eUseSegPar=UseSegPar;
    vertexrec->eQualityMode=QualityMode;

    vertex = dproc->GetVertexFromTree(*vertexrec,vertexfilename,vID);
  }
  else{
    vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");
    vertex = (EdbVertex*) vertexrec->eVTX->At(vID);
  }

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

void drawEDAvertices(bool newversion = true, TString vertexfilename= "vertextree_firstquarter_noendend.root"){
 using namespace VERTEX_PAR;
 TFile * inputfile = TFile::Open(vertexfilename.Data());
 TTree *vtxtree = (TTree*) inputfile->Get("vtx");
 EdbVertexRec *vertexrec;
 EdbDataProc *dproc = new EdbDataProc();


 TObjArray *drawnvertices = new TObjArray(1000000);
 TObjArray *drawntracks = new TObjArray(1000000);

// EdbTrackP *specialtrack = new EdbTrackP();

 EdbPVRec *ali = new EdbPVRec();
 EdbScanCond *scancond = new EdbScanCond();
 ali->SetScanCond(scancond);

 const int nvertices = 10000;
 cout<<"Reading number of vertices: "<<nvertices<<endl;
 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = i;
  
  EdbVertex *vertex = 0;

  if (newversion){
    vertexrec = new EdbVertexRec();
    vertexrec->SetPVRec(ali);
    vertexrec->eDZmax=DZmax;
    vertexrec->eProbMin=ProbMinV;
    vertexrec->eImpMax=ImpMax;
    vertexrec->eUseMom=UseMom;
    vertexrec->eUseSegPar=UseSegPar;
    vertexrec->eQualityMode=QualityMode;
    vertex = dproc->GetVertexFromTree(*vertexrec,vertexfilename,vID);
  } 
  else{ 
    vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");
    vertex = (EdbVertex*) vertexrec->eVTX->At(vID);
  }
  if (vertex->N()< 4) continue;
  drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
  for (int itrk = 0; itrk < vertex->N(); itrk++){
     EdbTrackP* track =  vertex->GetTrack(itrk);          
 //    for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[i]); // to color them differently
     drawntracks->Add(track);
//     specialtrack = track;
  }
 }

 ali->eTracks = drawntracks;
 ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
 eda->Run();
}
