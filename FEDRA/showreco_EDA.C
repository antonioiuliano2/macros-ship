//first script to launch automatic shower reconstruction from a list of tracks. Using EDA as interface to libShower (21 January 2020)

 EdbPVRec * dummyrec = new EdbPVRec();
 EdbEDA * eda = new EdbEDA(dummyrec);
 EdbEDAShowerTab *showtab = new EdbEDAShowerTab();
 TCanvas *cshow = new TCanvas(); 
 TGraph *memorygraph = new TGraph();
 MemInfo_t memorymonitor;
 ProcInfo_t procmonitor; 

void showreco_EDA(int trackID, EdbPVRec ** map_pvrec,fstream &logfile,TTree* tracks, int ipoint);
void showloop(){
  cshow->Divide(1,2);
  fstream inputtrackslist;
  inputtrackslist.open("inputtracklist.dat",fstream::in);
  gEDBDEBUGLEVEL = 0; //suppress FEDRA printout
  const int npvrecs = 30;
  EdbPVRec *map_pvrec[npvrecs];
  int ipvr = 0;
  cout<<"Getting map of PVRs from files"<<endl;
  for (int ix = 6; ix < 12; ix++){
      for (int iy = 0; iy < 5;iy++){
          int xstart = ix * 10000;
          int ystart = iy * 10000;
          TFile *myfile = TFile::Open(Form("pvrecs/pvrec_%d_%d.root",xstart,ystart));
          map_pvrec[ipvr] = (EdbPVRec*) myfile->Get("EdbPVRec");
          cout<<ipvr<<" "<<ix<<" "<<iy<<endl;
          ipvr++;//increasing counter
      }
  }
 cout<<"Finishing getting map"<<endl;
 
//tracks tree file
 TFile *trackfile = TFile::Open("linked_tracks.root");
 TTree *tracks = (TTree*)trackfile->Get("tracks");

 //const int ntracks_show = 10;
 fstream logfile;
 logfile.open("showrecotest.dat",fstream::out);
 //int inputtracks[ntracks_show] = {53650,16681,16683,8523,8520,55859,8505,86529,63838,129395};
// for (int itrk_show = 0; itrk_show<ntracks_show; itrk_show++){
 int ipoint = 0;
 while(inputtrackslist.good()){
  int trackID;
  inputtrackslist>>trackID;
  showreco_EDA(trackID,map_pvrec,logfile,tracks,ipoint);
  eda->Reset();
  ipoint++;
//  cout<<test<<endl;
 }
 TCanvas *cmemory = new TCanvas();
 memorygraph->Draw("AP*");
}

void showreco_EDA(int trackID, EdbPVRec ** map_pvrec,fstream &logfile,TTree* tracks, int ipoint){
    //array of EdbPVRec, opening files
 const int ntracks = 1;
 int tracklist = trackID; //87251
 
  //reading track file and setting branches
 int code = 0;
 int startplate = 0;

 Int_t   trid=0;
 Int_t   nseg=0;
 Int_t   npl=0;
 Int_t   n0=0;
 Float_t xv=0.;
 Float_t yv=0.;

 TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
 TClonesArray *segf = new TClonesArray("EdbSegP", 60);
 EdbSegP *trk=0;
 EdbSegP *s1=0;
 EdbSegP *s1f=0;

 tracks->SetBranchAddress("trid", &trid);
 tracks->SetBranchAddress("nseg", &nseg);
 tracks->SetBranchAddress("npl", &npl);
 tracks->SetBranchAddress("n0", &n0);
 tracks->SetBranchAddress("xv", &xv);
 tracks->SetBranchAddress("yv", &yv);
 tracks->SetBranchAddress("sf", &segf);
 tracks->SetBranchAddress("s",  &seg);
 tracks->SetBranchAddress("t.", &trk);

 //Getting Entry

 tracks->BuildIndex("trid");
 EdbPattern *pat = 0;
  //if (!(ali->GetTrack(trackID))) continue;
  seg->Clear();
  segf->Clear();
  tracks->GetEntryWithIndex(trackID);
  

  EdbTrackP *tr1 = new EdbTrackP();
  ((EdbSegP*)tr1)->Copy(*trk);
  tr1->SetM(0.139);                 //TODO  
    for(int i=0; i<nseg; i++) {
      s1 = (EdbSegP*)(seg->At(i));
      if(i == 0){
      int xcode = (s1->X()-60000)/10000;
      int ycode = s1->Y()/10000;
      code = 5*xcode + ycode;
      cout<<"POSITION "<<s1->X()<<" "<<s1->Y()<<endl;
      cout<<" TEST CODING "<<xcode<< " "<<ycode<<" "<<code<<endl;
      if (code > 30) code = 30;

      }
      pat = map_pvrec[code]->GetPattern(s1->PID());
      if (!pat){
        cout<<"Creating new pattern"<<endl;
        pat = new EdbPattern( 0., 0., s1->Z());
        pat->SetID(s1->PID());
        pat->SetScanID(s1->ScanID());
        map_pvrec[code]->AddPatternAt(pat,s1->PID());
      }

      s1f = (EdbSegP*)(segf->At(i));
      s1->SetDZ(300);
      s1f->SetDZ(300);
      s1->SetFlag(kMagenta); //EDA colours tracks according to flag

      EdbSegP * news1 = new EdbSegP(*((EdbSegP*)(s1)));
      EdbSegP * news1f = new EdbSegP(*((EdbSegP*)(s1f)));
      tr1->AddSegment( news1  );
      tr1->AddSegmentF( news1f);
    } //end on loop on segments associated to the track
  tr1->SetSegmentsTrack(tr1->ID());
  tr1->SetCounters();    

  //adding track to drawing set and building EDA
  map_pvrec[code]->AddTrack(tr1);

  startplate = tr1->GetSegmentFirst()->Plate(); 
// commenting eda part to test memory leakage
 eda->SetPVR(map_pvrec[code]);
 eda->GetTrackSet("TS")->AddTracksPVR(map_pvrec[code]);
 //eda = new EdbEDA(ali);
// EdbEDAUtil::FillTracksFromPatterns(ali);
 eda->GetTrackSet("TS")->SetColorMode(kCOLOR_BY_PARTICLE);
 eda->Run();
// eda->SavePictures();
 TObjArray * segarray = new TObjArray(100);
 segarray->Add(eda->GetTrackSet("TS")->GetTrack(0)->GetSegmentFirst());
 eda->SetSelected(segarray);


// showtab->Reco();

 //monitoring meminfo
 gSystem->GetMemInfo(&memorymonitor);
 gSystem->GetProcInfo(&procmonitor);

 cout<<"MEMORY MONITOR: "<<memorymonitor.fMemUsed<<" PROC MONITOR "<<procmonitor.fMemResident<<endl;

 memorygraph->SetPoint(ipoint,ipoint,procmonitor.fMemResident);

 TFile* file=TFile::Open("Shower.root");
 TTree* ShowerTree = (TTree *)file->Get("treebranch");

 int showersize;
 float output15, output30;

 ShowerTree->SetBranchAddress("sizeb",&showersize);
  ShowerTree->SetBranchAddress("output15",&output15);
 ShowerTree->SetBranchAddress("output30",&output30);
 cshow->cd(1);
 ShowerTree->Draw("yb:zb","","*");
 cshow->cd(2);
 ShowerTree->Draw("xb:zb","","*");
 cshow->Print(Form("showpics/showpicture_%d.png",trackID));
 int nshowers = ShowerTree->GetEntries();
 //logfile<<"Shower Selector: Segment from plate"<<startplate<<endl;
 if (nshowers == 0) cout<<"No shower found!"<<endl;
 else{
   //logfile<<"Startplate\tShower size\toutput15\toutput30"<<endl;
   for (int ishower = 0; ishower < nshowers; ishower++){
     ShowerTree->GetEntry(ishower);
     logfile<<startplate<<"\t"<<showersize<<"\t"<<output15<<"\t"<<output30<<endl;
   }



 }
 file->Close();

 //finishing showreco for this track, removing tracks from PVRec:
//  map_pvrec[code]->eTracks = 0;

//  delete segarray;
//  delete showtab;

 map_pvrec[code]->ResetTracks();
  delete seg;
  delete segf;
  delete tr1;
 //delete eda;
 //delete ali;
}