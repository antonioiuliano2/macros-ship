void drawEDAvertex(TString vertexfilename= "vertices_MC_modified.root"){
 const int nvertices = 2;
 int vertexlist[nvertices] = {25874,25879};

 TFile * inputfile = TFile::Open(vertexfilename.Data());
 EdbVertexRec *vertexrec = (EdbVertexRec*) inputfile->Get("EdbVertexRec");


 TObjArray *drawnvertices = new TObjArray(100);
 TObjArray *drawntracks = new TObjArray(100);

 for (int vID: vertexlist){ //range for loop, C++11
  EdbVertex *vertex = (EdbVertex*) vertexrec->eVTX->At(vID);

  drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
  for (int itrk = 0; itrk < vertex->N(); itrk++){
     EdbTrackP* track =  vertex->GetTrack(itrk);
     drawntracks->Add(track);
  }
 }
 EdbPVRec *ali = new EdbPVRec();

 ali->eTracks = drawntracks;
 ali->eVTX = drawnvertices;
 EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
 eda->Run();
}
