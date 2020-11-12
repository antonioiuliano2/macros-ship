//A new code to add protons to reconstructed vertices, after discussion with Giuliana Galati (GG) on 29 October 2020
void set_segments_dz(float dz);

EdbPVRec *ali = new EdbPVRec();

float CalcIP(EdbTrackP *tr, TVector3 V){
        //transverse distance (IP) from track to vertex (n.d.r. tranvserse with respect to beam z direction),
        //taken from nearest upstream segment
        float imp = 10000.;
        for (int iseg = 0; iseg < tr->N(); iseg++){
           EdbSegP *seg = tr->GetSegmentF(iseg);
           if (seg->Z() <= V(2)){           
            float dz = V(2) - seg->Z();
            float ipx = seg->TX() * dz + seg->X() - V(0);
            float ipy = seg->TY() * dz + seg->Y() - V(1);
            imp = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));           
           }
         } 
        return imp;

}

bool SplitTrackAtZ( EdbTrackP &t, EdbTrackP &t1, float vz )
{
  // Code copied from EdbTrackFitter::SplitTrack to split a track in 2. 
  // Only here I set the split point at z of vertex, instead of in a certain segment number
  for(int i=t.N()-1; i>=0; i--) {
    EdbSegP *seg = t.GetSegment(i);
    if (seg->Z() > vz){
     t1.AddSegment(   t.GetSegment(i) );
     t1.AddSegmentF(   t.GetSegmentF(i) );
     t.RemoveSegment( t.GetSegment(i) );
    }
  }
  t.SetCounters();
  //t1.FitTrackKFS(true, X0, 0);
  t1.SetCounters();
  t1.SetM(t.M());
  t1.SetP(t.P());
  //t1.FitTrackKFS(true, X0, 0);

  return true;
}

void addprotons(TString vertexfilename= "vertextree_firstquarter.root"){
    //selection which good vertices use as sample
    float dodisplay = false;
    TCut vertexselection = TCut("n>=10");
    const float maximp = 10; //cut on IP to accept the track (to be set  from MC)
    int nfound = 0;

    TObjArray *drawnvertices = new TObjArray(1000000);
    TObjArray *drawntracks = new TObjArray(1000000);

    EdbDataProc *dproc = new EdbDataProc();

    TH1D *hip = new TH1D("hip","IP of vertex with closest passing track;IP[#mum]",20,0,20);
    TH1I *hMCcheck = new TH1I("hMCcheck","Did we find the right particle? MC Check",2,0,2);
    TProfile *hMCcheck_ip = new TProfile("hMCcheck_ip","Did we find the right particle? MC Check vs IP",20,0,20,0,2);
    //track subset, 3 micron around it in xy distance
    TCut tracksel("nseg>3&&npl>=5&&t.Theta()<0.05");
    dproc->InitVolume(100, tracksel);
    ali = dproc->PVR();
    ali->FillCell(30,30,0.009,0.009);

    int nfoundtrue = 0;
    //adding vertexrec
    EdbVertexRec *vrec = new EdbVertexRec();
    vrec->SetPVRec(ali);
    //preparing the grid, to select indexes only from grid
    const int ncellsX = 12;
    const int ncellsY = 10;
    const float xmin = 0;
    const float xmax = 125000;
    const float ymin = 0;
    const float ymax = 100000;
    EdbH2 xygrid(ncellsX, xmin, xmax, ncellsY, ymin, ymax);
    vector<int> gridtracks[ncellsX][ncellsY];

    //filling tracks in grid;
    for (int itrk = 0; itrk < ali->eTracks->GetEntries(); itrk++){
      EdbTrackP * tr = (EdbTrackP*) ali->eTracks->At(itrk);
      int celltrx = xygrid.IX(tr->GetSegmentFirst()->X());
      int celltry = xygrid.IY(tr->GetSegmentFirst()->Y());
      gridtracks[celltrx][celltry].push_back(itrk);
      xygrid.Fill(tr->GetSegmentFirst()->X(), tr->GetSegmentFirst()->Y());
    }
    
    cout<<"Track grid was filled with size "<<ncellsX<<" times "<<ncellsY<<endl;
    //getting vertex
    dproc->ReadVertexTree(*ali,vertexfilename,vertexselection);
    const int nvertices = ali->eVTX->GetEntries();
    set_segments_dz(300.);

    //we want to know the most frequent MCEvent within vertex
    map<int,int> frequencyEvent;

    cout<<"Start loop on vertices"<<endl;
    fstream outputlist("10tracksvertices_addedprotons.log",fstream::out);
    for (int ivertex = 0; ivertex < nvertices; ivertex++){
     int foundtrueproton = 0;
     frequencyEvent.clear();

     EdbVertex *vertex = (EdbVertex*) ali->eVTX->At(ivertex);
     float vx = vertex->X();
     float vy = vertex->Y();
     float vz = vertex->Z();

     int vID = vertex->ID();

     //adding vertex and its tracks to be drawn
     drawnvertices->Add(vertex); // assuming the array is filled with EdbVertex.
     float minimp = maximp;
     EdbTrackP *protoncandidate;
     //finding most frequent MCEvent in this vertex
     for (int itrk = 0; itrk < vertex->N(); itrk++){
  
      EdbTrackP* track =  vertex->GetTrack(itrk);
      int eventtrack = track->GetSegmentFirst()->MCEvt();
      frequencyEvent[eventtrack]++;

      for (int iseg=0; iseg<track->N();iseg++) track->GetSegment(iseg)->SetFlag(kRed); //standard color for initial particles
     //for (int iseg = 0; iseg < track->N(); iseg++) track->GetSegment(iseg)->SetFlag(vertexcolors[i]); // to color them differently
      drawntracks->Add(track);
      //specialtrack = track;
     }
     //finding most common event
     map<int,int>::iterator it; 
     int ntracks_event = 0;
     int mostfrequentevent = -1;
     for (it = frequencyEvent.begin(); it!=frequencyEvent.end();it++){
       if(it->second > ntracks_event){
        ntracks_event = it->second;
        mostfrequentevent = it->first;
       }
     }
     //which tracks do we need to loop into?
     int cellvx = xygrid.IX(vx);
     int cellvy = xygrid.IY(vy);
     //start main loop on external tracks
     for (int itrk = 0; itrk < gridtracks[cellvx][cellvy].size(); itrk++){
      int trackindex = gridtracks[cellvx][cellvy].at(itrk);
      EdbTrackP * tr = (EdbTrackP*) ali->eTracks->At(trackindex);
        if (tr->GetSegmentFirst()->Z() < vertex->Z()){ //proton candidate should start before vz
         
         //transverse IP to vertex
         float imp =  CalcIP(tr, TVector3(vx,vy,vz));
        
        
         if(imp < minimp){ //found candidate, storing track and new min IP
          protoncandidate = tr;
          minimp = imp;          
          } //imp condition
         } //first trk condition

        }//loop on trks
        hip->Fill(minimp);
        if (minimp < maximp){ 
            for (int iseg=0; iseg<protoncandidate->N();iseg++) protoncandidate->GetSegment(iseg)->SetFlag(kGreen); //i label the candidates differently
            
            

            if (protoncandidate->GetSegmentLast()->Z() > vz){ //track ending after vertex, splitting into two
             EdbTrackP * protondaughter = new EdbTrackP();
             SplitTrackAtZ(*protoncandidate, *protondaughter, vz);

             drawntracks->Add(protondaughter);
             EdbVTA *vta1 = new EdbVTA(protondaughter,vertex);
             vta1->SetFlag(2);
             vertex->AddVTA(vta1);
             vta1->SetZpos(1);
             protoncandidate->AddVTA(vta1);
            }

            drawntracks->Add(protoncandidate);
            EdbVTA *vta = new EdbVTA(protoncandidate,vertex);
            vta->SetFlag(2);
            vertex->AddVTA(vta);
            vta->SetZpos(0);
            protoncandidate->AddVTA(vta);

            nfound++;                  
            if (protoncandidate->MCEvt() == mostfrequentevent && protoncandidate->GetSegmentFirst()->MCTrack()==0) foundtrueproton = 1;
          
            outputlist << vID << " " << protoncandidate->Track() << " "<<endl;
        }
      hMCcheck_ip->Fill(minimp, foundtrueproton);   
      hMCcheck->Fill(foundtrueproton); 

      if (foundtrueproton) nfoundtrue++;
    } //loop on vertices    
    cout<<"Found candidates protons: "<<nfound<<" true MC"<<nfoundtrue<<endl;
    
    if(dodisplay){
     ali->eTracks = drawntracks;
     ali->eVTX = drawnvertices;
     EdbEDA * eda = new EdbEDA(ali); // init DataSet but doesn't read linked_track.root
     eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
     eda->Run();
    }

    TCanvas *c = new TCanvas();
    hip->Draw();

    TCanvas *cMC = new TCanvas();
    hMCcheck->Draw();
    
    TCanvas *cMC2D = new TCanvas();
    hMCcheck_ip->Draw();

    TCanvas *ctrackgrid = new TCanvas();
    xygrid.DrawH2("hxygrid","Grid with tracks, from first segment position;x[#mu m];y[#mu m]")->Draw("COLZ");

    outputlist.close();
}

//---------------------------------------------------------------------
void set_segments_dz(float dz)
{
  int np = ali->Npatterns();
  for(int i=0; i<np; i++) {
    EdbPattern *p = ali->GetPattern(i);
    int ns = p->N();
    for(int j=0; j<ns; j++) p->GetSegment(j)->SetDZ(dz);
  }
}