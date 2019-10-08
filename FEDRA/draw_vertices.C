/*script for vertex file analysis
By default does drawing of the vertices and of the plots in a unique loop, along with a production of a dump file
root -l
.L analysis_vertices()
analysis_vertices(ivertex) to draw vertex with ID ivertex
*/
TString inputfilename = "vertices_MC_modified.root";

//analyze and draw vertexes with multiplicity>=trmin, and aperture >= amin

void draw_selected_vertices(int trmin=2, float amin=0.01){
 const int nvertices = 2;
 int vertexlist[nvertices] = {25874,25879};
 
 TFile *inputfile = TFile::Open(inputfilename.Data()); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;
 EdbVertexRec *gEVR = (EdbVertexRec*) inputfile->Get("EdbVertexRec");

 //arrays for vertices drawing
 TObjArray *varr = new TObjArray();
 TObjArray *tarr = new TObjArray();

 EdbVertex *v=0;
 EdbTrackP *t1=0;

 for (int vID: vertexlist){
    v = (EdbVertex *)(gEVR->eVTX->At(vID));
    int ntracks = v->N();
    cout<<"Flag: "<<v->Flag()<<" ntracks: "<<v->N()<<" maxaperture: "<<v->MaxAperture()<<endl;
    //cuts on vertices    
    if(v->Flag()<0)         continue;
    if( ntracks<trmin) continue;
    if( v->MaxAperture()<amin )  continue;

    varr->Add(v); //adding vertex to array for drawing

    //LOOP ON TRACKS
    for(int j=0; j<ntracks; j++){

     //acceding to the track
     t1 = v->GetTrack(j);
  
     tarr->Add(t1); //adding track to array for drawing        
   } //end of tracks loop
 }//end of vertex loop
  const char *dsname="display-v";
  EdbDisplay   *ds=0;
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

void draw_vertices(int vid = -1,int trmin=2, float amin=0.01){

  TFile *inputfile = TFile::Open(inputfilename.Data()); 
  if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;

  EdbVertexRec *gEVR = (EdbVertexRec*) inputfile->Get("EdbVertexRec");

  //arrays for vertices drawing
  TObjArray *varr = new TObjArray();
  TObjArray *tarr = new TObjArray();

  EdbVertex *v=0;
  EdbTrackP *t1=0;

  int nv = gEVR->Nvtx();
  printf("Number of vertices =%d\n",nv);
  if(nv<1) return;

  //*********************LOOP ON VERTICES****************************
  for(int i=0; i<nv; i++) {
    if ((vid != -1) && (i != vid)) continue;

    v = (EdbVertex *)(gEVR->eVTX->At(i));
    int ntracks = v->N();
    cout<<"Flag: "<<v->Flag()<<" ntracks: "<<v->N()<<" maxaperture: "<<v->MaxAperture()<<endl;
    //cuts on vertices    
    if(v->Flag()<0)         continue;
    if( ntracks<trmin) continue;
    if( v->MaxAperture()<amin )  continue;

    varr->Add(v); //adding vertex to array for drawing

    //LOOP ON TRACKS
    for(int j=0; j<ntracks; j++){

     //acceding to the track
     t1 = v->GetTrack(j);
  
     tarr->Add(t1); //adding track to array for drawing        
   } //end of tracks loop
  } //end of vertex loop

  const char *dsname="display-v";
  EdbDisplay   *ds=0;
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
