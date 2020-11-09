//display for vertices, patterns are not drawn

void drawEDAvertex(bool newversion = true, TString vertexfilename= "vertextree.root"){
 const int nvertices = 1;
 int vertexlist[nvertices] = {110839};
 int vertexcolors[nvertices] = {kRed};
 EdbDataProc *dproc = new EdbDataProc();

 TFile * inputfile = TFile::Open(vertexfilename.Data());
 EdbVertexRec *vertexrec; 

 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);

// EdbTrackP *specialtrack = new EdbTrackP();

 EdbPVRec *ali = new EdbPVRec();

 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = vertexlist[i];
  
  EdbVertex *vertex = 0;

  if (newversion) vertex = dproc->GetVertexFromTree(*ali,vertexfilename,vID);
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
 TFile * inputfile = TFile::Open(vertexfilename.Data());
 TTree *vtxtree = (TTree*) inputfile->Get("vtx");
 EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");
 EdbDataProc *dproc = new EdbDataProc();


 TObjArray *drawnvertices = new TObjArray(1000000);
 TObjArray *drawntracks = new TObjArray(1000000);

// EdbTrackP *specialtrack = new EdbTrackP();

 EdbPVRec *ali = new EdbPVRec();

 const int nvertices = 10000;
 cout<<"Reading number of vertices: "<<nvertices<<endl;
 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = i;
  
  EdbVertex *vertex = 0;

  if (newversion) vertex = dproc->GetVertexFromTree(*ali,vertexfilename,vID);
  else vertex = (EdbVertex*) vertexrec->eVTX->At(vID);

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
