void drawEDAvertex(int vID = 0){
 const int nvertices = 1;
 int vertexlist[nvertices] = {3781};

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
 EdbPVRec *pvr = new EdbPVRec();
 pvr->eTracks = drawntracks;
 pvr->eVTX = drawnvertices;

 EdbEDA * eda = new EdbEDA(pvr); // init DataSet but doesn't read linked_track.root
 eda->Run();

}
