//display for vertices, patterns are not drawn
void drawEDAvertex(TString vertexfilename= "vertices_MC_modified.root"){
 const int nvertices = 3;
 int vertexlist[nvertices] = {25578,25972,24660};
 int vertexcolors[nvertices] = {kRed,kGreen,kBlue};

 TFile * inputfile = TFile::Open(vertexfilename.Data());
 EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");


 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);

// EdbTrackP *specialtrack = new EdbTrackP();

 for (int i = 0; i < nvertices; i++){ //range for loop, C++11
  int vID = vertexlist[i];
  EdbVertex *vertex = (EdbVertex*) vertexrec->eVTX->At(vID);

  drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
  for (int itrk = 0; itrk < vertex->N(); itrk++){
     EdbTrackP* track =  vertex->GetTrack(itrk);
     for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[i]); // to color them differently
     drawntracks->Add(track);
//     specialtrack = track;
  }
 }
 EdbPVRec *ali = new EdbPVRec();

 ali->eTracks = drawntracks;
 ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
/* eda->GetTrackSet("TS")->RemoveTrack(specialtrack); //charm colored differently
 eda->GetTrackSet("SB")->SetTrackAttribute(4);
 eda->GetTrackSet("SB")->AddTrack(specialtrack);*/
 eda->Run();
}
