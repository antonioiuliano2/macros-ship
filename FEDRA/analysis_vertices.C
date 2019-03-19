/*script for vertex file analysis
By default does drawing of the vertices and of the plots in a unique loop, along with a production of a dump file
root -l
.L analysis_vertices()
analysis_vertices(ivertex) to draw vertex with ID ivertex
*/
TString run = "CH1-R6";
TString path = "/ship/CHARM2018/" + run +"/b000001/"; 

TString inputfilename = path + "vertices.root";

//analyze and draw vertexes with multiplicity>=trmin, and aperture >= amin
void analysis_vertices(int vid = -1,int trmin=4, float amin=0.01){

  TFile *inputfile = TFile::Open(inputfilename.Data()); 
  if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;

  EdbVertexRec *gEVR = (EdbVertexRec*) inputfile->Get("EdbVertexRec");

  //histograms for checks of kinematics

  TH2D *hvprobn = new TH2D("hvprobn", "Vertex Quality vs molteplicity",50,0,50,30,0,0.3);
  TH2D *hvxy = new TH2D("hvxy","vertex position", 120, 0, 120000, 100, 0, 100000);
  TH2D *hfill_n = new TH2D("hfill_n", "Fill factor in quadratic mean for tracks associated to vertex",50,0,50,20,0,1);
  TH2D *hdeltathetaeff = new TH2D("hdeltathetaeff", "Efficiency vs maximum difference in theta angle",100,0,0.1, 100, 0, 1.0);
 
  TH1D *hvz = new TH1D("hvz","vertex position z", 400, -40000, 0);
  TH1D *hvn = new TH1D("hvn","vertex molteplicity", 50, 0, 50);
  TH1D *hip = new TH1D("hip","impact parameter of tracks with respect to vertex",100,0,100);
  TH1D *htheta = new TH1D("htheta", "Theta angle tracks", 3500, 0, 3.5);

  TH1D *hdeltatheta = new TH1D("hdeltatheta", "Difference in teta angle between volume track angle and segment angle", 1000, 0, 1.0);

  TH1D *hrmslongaverage = new TH1D("hrmslongaverage","Longitudinal", 100, 0, 0.1);
  TH1D *hrmstransaverage = new TH1D("hrmstransaverage","Transverse",100, 0, 0.1);
  TH1D *hrmslongmany = new TH1D("hrmslongmany","Longitudinal", 100, 0, 0.1);
  TH1D *hrmstransmany = new TH1D("hrmstransmany","Transverse",100, 0, 0.1);
  TH1D *hrmslongfew = new TH1D("hrmslongfew","Longitudinal", 100, 0, 0.1);
  TH1D *hrmstransfew = new TH1D("hrmstransfew","Transverse",100, 0, 0.1);
  TH1D *hrmsspaceaverage = new TH1D("hrmsspaceaverage","Space angle", 100, 0, 0.1);
  TH1D *hrmsspacemany = new TH1D("hrmsspacemany","Space angle", 100, 0, 0.1);
  TH1D *hrmsspacefew = new TH1D("hrmsspacefew","Space angle", 100, 0, 0.1);
  //arrays for vertices drawing
  TObjArray *varr = new TObjArray();
  TObjArray *tarr = new TObjArray();

  EdbVertex *v=0;
  EdbTrackP *t1=0;

  int nv = gEVR->Nvtx();
  printf("Number of vertices =%d\n",nv);
  Int_t ntracksstart, ntracksend;
  if(nv<1) return;

  Double_t track_fill = 0.;
  Double_t meanfill;

  //*********************LOOP ON VERTICES****************************
  for(int i=0; i<nv; i++) {
    if ((vid != -1) && (i != vid)) continue;
    ntracksstart = 0;
    ntracksend = 0;
    meanfill = 0.;

    v = (EdbVertex *)(gEVR->eVTX->At(i));
    int ntracks = v->N();
    //cout<<v->Flag()<<" "<<v->N()<<" "<<v->MaxAperture()<<" "<<v->VX()<<" "<<v->VY()<<endl;
    //cuts on vertices    
    if(v->Flag()<0)         continue;
    if( ntracks<trmin) continue;
    if( v->MaxAperture()<amin )  continue;
   /* for(int j=0; j<ntracks; j++){
        Int_t zpos = v->GetVTa(j)->Zpos();//ask if track is start or end: 0 track ends at vertex, 1 track start at vertex 
        if (zpos == 0) ntracksend++;
	if (zpos == 1) ntracksstart++;    
    }
    if (ntracksend != 1) continue;*/
    //end of cuts
    //filling histograms with vertices informations
    //cout<<"Probabilità del vertice: "<<v->V()->prob()<<endl;
    hvprobn->Fill(ntracks,v->Quality());
    hvxy->Fill(v->VX(), v->VY()); 
    hvz->Fill(v->VZ());
    hvn->Fill(ntracks);
    varr->Add(v); //adding vertex to array for drawing

    //LOOP ON TRACKS
    for(int j=0; j<ntracks; j++){

     //acceding to the track
     hip->Fill(v->Impact(j));
     t1 = v->GetTrack(j);

     //Defining the variable for studying the kink angle
     double rmsspace = 0.; 
     double rmstransverse = 0.;
     double rmslongitudinal = 0.;
     int NKinkAngleUsed = 0;
     double pippo = EdbEDAUtil::DTRMSTL1Kink(t1, &rmsspace, &rmstransverse, &rmslongitudinal, &NKinkAngleUsed); //prova se funziona il metodo
     //Filling the histograms

     if (ntracks<5){
     hrmstransfew->Fill(rmstransverse);
     hrmslongfew->Fill(rmslongitudinal);
     hrmsspacefew->Fill(rmsspace);
     }
     else if (ntracks>10){
     hrmstransmany->Fill(rmstransverse);
     hrmslongmany->Fill(rmslongitudinal);
     hrmsspacemany->Fill(rmsspace);
     }
     else{
     hrmstransaverage->Fill(rmstransverse);
     hrmslongaverage->Fill(rmslongitudinal);
     hrmsspaceaverage->Fill(rmsspace);
     }
     tarr->Add(t1); //adding track to array for drawing     

     track_fill = (Double_t)(t1->N())/(t1->Npl());//computing the 'fill factor' for a single track
     meanfill += track_fill;     

     Double_t theta1 = t1->Theta();
     //filling histograms with tracks informations
     htheta->Fill(theta1);
     //LOOP on SEGMENTS
     Double_t maxdeltatheta = 0.;
     for (int k=0; k<t1->N();k++){
	EdbSegP* seg = t1->GetSegment(k);
	Double_t thetaseg = seg->Theta();
        if (TMath::Abs(thetaseg-theta1)>maxdeltatheta) maxdeltatheta = TMath::Abs(thetaseg-theta1);
	hdeltatheta->Fill(TMath::Abs(thetaseg-theta1));
     }//end of segment loop
   hdeltathetaeff->Fill(maxdeltatheta, track_fill);
   } //end of tracks loop
  meanfill = meanfill/ntracks;
  hfill_n->Fill(ntracks,meanfill);
  } //end of vertex loop

 if (vid == -1){ //useless to draw plots for only one vertex
  //**********Drawing canvas************************************
  TCanvas * cxy = new TCanvas();
  hvxy->GetXaxis()->SetTitle("x[#mu m]");
  hvxy->GetYaxis()->SetTitle("y[#mu m]");
  hvxy->Draw();
 
  TCanvas * cz = new TCanvas();
  hvz->GetXaxis()->SetTitle("z[#mu m]");
  hvz->Draw();
  
  TCanvas *crmscomponents = new TCanvas();
  crmscomponents->Divide(1,2);
  crmscomponents->cd(1);
  hrmslongfew->Scale(1./hrmslongfew->GetEntries());
  hrmslongmany->Scale(1./hrmslongmany->GetEntries());
  hrmslongmany->SetLineColor(kRed);
  hrmslongaverage->Scale(1./hrmslongaverage->GetEntries());
  hrmslongaverage->SetLineColor(kYellow);
  hrmslongaverage->GetXaxis()->SetTitle("rad");
  hrmslongmany->Draw();
  hrmslongaverage->Draw("SAMES");
  hrmslongfew->Draw("SAMES");
  TLegend* legend0 = new TLegend(0.1,0.7,0.48,0.9);
  legend0->AddEntry("hrmslongfew", "Vertices with less than 5 tracks");
  legend0->AddEntry("hrmslongmany","Vertices with more than 10 tracks");
  legend0->AddEntry("hrmslongaverage","Vertices with tracks between 5 and 10");
  legend0->Draw("SAME");
  crmscomponents->cd(2);
  hrmstransfew->Scale(1./hrmstransfew->GetEntries());
  hrmstransmany->Scale(1./hrmstransmany->GetEntries());
  hrmstransmany->SetLineColor(kRed);
  hrmstransaverage->Scale(1./hrmstransaverage->GetEntries());
  hrmstransaverage->SetLineColor(kYellow);
  hrmstransaverage->GetXaxis()->SetTitle("rad");
  hrmstransmany->Draw();
  hrmstransfew->Draw("SAMES");
  hrmstransaverage->Draw("SAMES");
  TLegend* legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry("hrmstransfew", "Vertices with less than 5 tracks");
  legend1->AddEntry("hrmstransmany","Vertices with more than 10 tracks");
  legend1->AddEntry("hrmstransaverage","Vertices with tracks between 5 and 10");
  legend1->Draw("SAME");

  TCanvas *crmsspace = new TCanvas();
  hrmsspacefew->Scale(1./hrmsspacefew->GetEntries());
  hrmsspacemany->Scale(1./hrmsspacemany->GetEntries());
  hrmsspacemany->SetLineColor(kRed);
  hrmsspaceaverage->Scale(1./hrmsspaceaverage->GetEntries());
  hrmsspaceaverage->SetLineColor(kYellow);
  hrmsspaceaverage->GetXaxis()->SetTitle("rad");
  hrmsspacemany->Draw();
  hrmsspaceaverage->Draw("SAMES");
  hrmsspacefew->Draw("SAMES");
  TLegend* legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry("hrmsspacefew", "Vertices with less than 5 tracks");
  legend2->AddEntry("hrmsspacemany","Vertices with more than 10 tracks");
  legend2->AddEntry("hrmsspaceaverage","Vertices with tracks between 5 and 10");
  legend2->Draw("SAME");

  TCanvas * cdeltathetaeff = new TCanvas();
  hdeltathetaeff->GetXaxis()->SetTitle("deltatheta[rad]");
  hdeltathetaeff->GetYaxis()->SetTitle("Efficiency");
  hdeltathetaeff->Draw("COLZ"); 

  TCanvas * cmol = new TCanvas();
  hvn->Draw();

  TCanvas * cip = new TCanvas(); 
  hip->GetXaxis()->SetTitle("IP[#mu m]");
  hip->Draw();

  TCanvas * ctheta = new TCanvas();
  htheta->GetXaxis()->SetTitle("theta[rad]");
  htheta->Draw();

  TCanvas * cdeltatheta = new TCanvas();
  hdeltatheta->GetXaxis()->SetTitle("deltatheta[rad]");
  hdeltatheta->Draw();

  TCanvas *cprob = new TCanvas();
  hvprobn->Draw("COLZ");
  hvprobn->GetXaxis()->SetTitle("NTracks");
  hvprobn->GetYaxis()->SetTitle("Vertex quality parameter [#mum^{-4}]");

  //Studying fill factor
  TCanvas *ctrackfill = new TCanvas();
  hfill_n->Draw("COLZ");  

/*//try with EDA
  eda = new EdbEDA("lnk.def", -1);
  eda->GetTrackSet("TS")->AddTracks(tarr);
  eda->GetVertexSet()->AddVertices(varr);
  eda->Run();
*/
  //vertex drawing
  gStyle->SetPalette(1);
 }

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


void dump_vertices(int trmin=4, float amin=0.01){

 TFile *inputfile = TFile::Open(inputfilename.Data()); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;
 EdbVertexRec *gEVR = (EdbVertexRec*) inputfile->Get("EdbVertexRec");

 TString dumpfilename = path + "vertexlist.txt";
 TString setfilename = path + "b000001.0.0.0.set.root";

 TFile *setfile = TFile::Open(setfilename.Data());
 EdbScanSet *scanset = (EdbScanSet*) setfile->Get("set");

 EdbVertex * vertex = new EdbVertex;
 int nv = gEVR->Nvtx();

 fstream outputfile(dumpfilename.Data(),fstream::out);

 Double_t Vx, Vy, Vz; //vertex position in global coordinates
 Double_t Vx_local, Vy_local; //vertex position in local coordinates
 Double_t Vertex_DZ;
 Double_t trackstartx_local,trackstarty_local,trackstartz;
 EdbSegP *firstsegment; //first segment of the track

 Int_t ID_emu_downstream; //number of the emulsion downstream
 Int_t ntracksstart, ntracksend;
 Int_t ntracksdownstream; //number of tracks associated to the vertex and found in downstream plate
 EdbAffine2D *aff; //affine transformation from local to global;
 EdbAffine2D *invertaff;

 outputfile<<"Format for track information: TrackID, (startx, starty), (TX,TY), track->N()"<<endl;
  for(int i=0;i<nv;i++){//accessing vertex list
    //COUNTER RESET
    ntracksstart = 0;
    ntracksend = 0;
    ID_emu_downstream = -1;    
    ntracksdownstream = 0;    
    Vertex_DZ = 0.;        
    ////////////////

    vertex = (EdbVertex *)(gEVR->eVTX->At(i));
    int ntracks = vertex->N();
    Vx = vertex->VX();
    Vy = vertex->VY();
    Vz = vertex->VZ();
    //*****************CUTS ON VERTEX
    //cuts on vertices    
    if(vertex->Flag()<0)         continue;
    if( ntracks<trmin) continue;
    if( vertex->MaxAperture()<amin )  continue;
/*    for(int j=0; j<ntracks; j++){
        Int_t zpos = vertex->GetVTa(j)->Zpos();//ask if track is start or end: 0 track ends at vertex, 1 track start at vertex 
        if (zpos == 0) ntracksend++;
	if (zpos == 1) ntracksstart++;    
    }    
    if (ntracksend != 1) continue;
    ntracksstart = 0;
    ntracksend = 0;*/
    //*****************END OF CUTS
    for(int j=0; j<ntracks; j++){ //loop on associated tracks
		Int_t zpos = vertex->GetVTa(j)->Zpos();//ask if track is start or end: 0 track ends at vertex, 1 track start at vertex 
                if (zpos == 0) ntracksend++;
		if (zpos == 1) ntracksstart++;       
      
                EdbTrackP *t1 = vertex->GetTrack(j);
		firstsegment = t1->GetSegmentFirst();
                trackstartz = firstsegment->Z();
		//I need to get the ID of the first plate downstream of the vertex
		if ((trackstartz - Vz) > 0 && (trackstartz - Vz) < 2000){ 
			ntracksdownstream++;
			ID_emu_downstream = firstsegment->Plate();
			Vertex_DZ = trackstartz - Vz;
			
			//getting start of track in local coordinates
			aff = scanset->Brick()->GetPlate(29-ID_emu_downstream)->GetAffineXY();    
 			invaff = (EdbAffine2D*) aff->Clone();
     			invaff->Invert();
     			trackstartx_local = invaff->Xtrans(firstsegment->X(),firstsegment->Y()); 
     			trackstarty_local = invaff->Ytrans(firstsegment->X(),firstsegment->Y());

			outputfile<<t1->ID()<<" "<<trackstartx_local<<" "<<trackstarty_local<<" "<<t1->TX()<<" "<<t1->TY()<<" "<<t1->N()<<endl;
			//outputfile<<"Track ID: "<<t1->ID()<<" with number of basetracks: "<<t1->N()<<endl;
			//outputfile<<"Start of tracks in local coordinates (x,y):"<<trackstartx_local<<" "<<trackstarty_local<<endl;
			//outputfile<<"Track angles (Tx,Ty):" <<t1->TX()<<" "<<t1->TY()<<endl;
			//outputfile<<"Number of grians in s[0]:" <<firstsegment->W()<<endl;
			}
//		cout<<trackstartz-Vz<<" "<<trackstartz<<" "<<Vz<<endl;
              }
  //  cout<<Vz<<" "<<ID_emu_downstream<<endl;
    if (ID_emu_downstream > -1){
     aff = scanset->Brick()->GetPlate(29-ID_emu_downstream)->GetAffineXY();
    // cout<<"Affine transformation (local to global)"<<endl;
    // aff->Print();
     invaff = (EdbAffine2D*) aff->Clone();
     invaff->Invert();
     Vx_local = invaff->Xtrans(Vx,Vy); //getting vertex position in local coordinates of the plate
     Vy_local = invaff->Ytrans(Vx,Vy);
    }
    outputfile<<"*********End of tracks: VERTEX INFORMATION************"<<endl;
    outputfile<<"ID Vertex: "<<i<<" Flag: "<<vertex->Flag()<<" N plate downstream: "<<ID_emu_downstream<<endl;
    outputfile<<"Distance from downstream plate"<<Vertex_DZ<<endl;
    if ((Vertex_DZ < 400) && (Vertex_DZ > 0)) outputfile<<"SMALL DZ"<<endl;
    if ((Vertex_DZ < 200) && (Vertex_DZ > 0) && (ID_emu_downstream == 26)) cout<<"Very SMALL"<<" "<<i<<endl;
    outputfile<<"N tracks upstream: "<<ntracksend<<" N tracks downstream: "<<ntracksstart<<endl;
    outputfile<<"N tracks in last film: "<<ntracksdownstream<<" N tracks total: "<<ntracks<<endl;
    outputfile<<"Global Coordinates: "<<Vx<<" "<<Vy<<" "<<Vz<<endl;
    outputfile<<"Local Coordinates: "<<Vx_local<<" "<<Vy_local<<" "<<Vz<<endl;
    outputfile<<endl;
 }
 outputfile.close();
}
