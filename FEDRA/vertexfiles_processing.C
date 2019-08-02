//routine to put in a simplest and clearest way possible all the needed steps after vertexing

//usual FEDRA functions to read from track tree

EdbDataProc  *dproc=0;
EdbPVRec     *gAli=0;

TString inputfilename = "vertices_MC.root"; //path to input vertex file

void trvol( const char *def, const char *rcut = "nseg>1" );
void init( const char *def, int iopt,  const char *rcut="1" );
void set_segments_dz(float dz=300.);
void savingremainingtracks(TObjArray* tracklist);
void modify_distribution_tree();

void trvol( const char *def, const char *rcut )
{
  // this function read volume tracks from linked_tracks.root

  init(def, 100 ,rcut);                      // read tracks (option 100)
  gAli->FillCell(30,30,0.009,0.009);
}

//---------------------------------------------------------------------
void init( const char *def, int iopt,  const char *rcut)
{
  if(!def)  dproc = new EdbDataProc();
  else      dproc = new EdbDataProc(def);

  dproc->InitVolume(iopt, rcut);  
  gAli = dproc->PVR();
  set_segments_dz(300.);
}

//---------------------------------------------------------------------
void set_segments_dz(float dz)
{
  int np = gAli->Npatterns();
  for(int i=0; i<np; i++) {
    EdbPattern *p = gAli->GetPattern(i);
    int ns = p->N();
    for(int j=0; j<ns; j++) p->GetSegment(j)->SetDZ(dz);
  }
}

void vertexfiles_processing(){ //script to fill vertex tree with various information and produce a tree with not associated tracks
 trvol(0,"nseg>1");
 TObjArray *tracklist = gAli->eTracks; 

 TFile *inputfile = TFile::Open(inputfilename.Data()); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;
 TTree *vertextree = (TTree*) inputfile->Get("vtx");
 
 const Int_t nvertices = vertextree->GetEntries();
 //defining variables for storing tree branches

 Int_t n, vID;
 Float_t vx, vy,vz;
 Float_t vtx_max_aperture, probability;
 //setting tree addresses
 vertextree->SetBranchAddress("vID",&vID);
 vertextree->SetBranchAddress("vx",&vx);
 vertextree->SetBranchAddress("vy",&vy);
 vertextree->SetBranchAddress("vz",&vz);
 vertextree->SetBranchAddress("probability",&probability);
 vertextree->SetBranchAddress("n",&n);

 //VARIABLE AND BRANCH preparations for output tree

 const Int_t maxdim = 1000;
 //big arrays for containers of track variables
 Int_t nseg[maxdim];
 //inserting MC truth information
 Int_t MCEventID[maxdim];
 Int_t MCTrackID[maxdim];
 Int_t MCMotherID[maxdim];
 double rmsthetatransverse[maxdim];
 double rmsthetalongitudinal[maxdim];
 //track variables
 Int_t TrackID[maxdim];
 Float_t TX[maxdim];
 Float_t TY[maxdim];
 Float_t trackfill[maxdim];
 Float_t impactparameter[maxdim];
 Int_t incoming[maxdim]; //coming from the vertex
 Int_t trk_num_holes[maxdim]; // number of missed segment in a track
 Int_t trk_max_gap[maxdim];   // max number of consecutive missed segment in a track

 TFile *outfile = new TFile("vertices_MC_modified.root","RECREATE");
 TTree *outvertextree = new TTree("vtx","Vertices");
 ////TTree *outqualitytree = qualitytree->CloneTree();

 //vertex variables
 outvertextree->Branch("vID",&vID,"vID/I");
 outvertextree->Branch("vx",&vx,"vx/F");
 outvertextree->Branch("vy",&vy,"vy/F");
 outvertextree->Branch("vz",&vz,"vz/F");
 outvertextree->Branch("vtx_max_aperture",&vtx_max_aperture,"vtx_max_aperture/F");
 outvertextree->Branch("probability",&probability,"probability/F");
 outvertextree->Branch("n",&n,"n/I");
 //track variables (they are array with the number of tracks as size)
 outvertextree->Branch("TrackID",&TrackID,"TrackID[n]/I");
 outvertextree->Branch("TX",&TX,"TX[n]/F");
 outvertextree->Branch("TY",&TY,"TY[n]/F");
 outvertextree->Branch("nseg",&nseg,"nseg[n]/I");
 //outvertextree->Branch("rmsthetalongitudinal",rmsthetalongitudinal,"rmsthetalongitudinal[n]/F");
 //outvertextree->Branch("rmsthetatransverse",rmsthetatransverse,"rmsthetatransverse[n]/F");
 outvertextree->Branch("trackfill",trackfill,"trackfill[n]/F");
 outvertextree->Branch("incoming",&incoming,"incoming[n]/I");
 outvertextree->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
 outvertextree->Branch("trk_num_holes",&trk_num_holes,"trk_num_holes[n]/I");
 outvertextree->Branch("trk_max_gap",&trk_max_gap,"trk_max_gap[n]/I");
 //inserting MCtrue information
 outvertextree->Branch("MCEventID", &MCEventID, "MCEventID[n]/I");
 outvertextree->Branch("MCTrackID",&MCTrackID,"MCTrackID[n]/I");
 outvertextree->Branch("MCMotherID",&MCMotherID,"MCMotherID[n]/I");
 
 //we need some Edb objects to add new information
 EdbVertexRec *vertexrec = (EdbVertexRec*)inputfile->Get("EdbVertexRec");
 EdbVertex *vertexobject = 0;
 EdbTrackP *track = 0;
 //Vertex    *vt = 0;  da chiedere ad Antonio
 //**************************************LOOP ON VERTICES*******************************
 for (int ivtx = 0; ivtx < nvertices; ivtx++){
  vertextree->GetEntry(ivtx);
  vertexobject = (EdbVertex *)(vertexrec->eVTX->At(vID));
  vtx_max_aperture = vertexobject->MaxAperture();
  
  //*************************************LOOP ON TRACKS**********************************
  for (int itrk = 0; itrk < n; itrk++){
 
   track = vertexobject->GetTrack(itrk);
   TrackID[itrk] = track->Track(); //NOT ID(), which is resetted! The eTrack attribute of EdbSegP now allows association to original root tree
   nseg[itrk] = track->N();
   //track fill factor definition
   if (track->Npl() <=2) trackfill[itrk] = 1.;
   else trackfill[itrk] = (Float_t)(track->N()-2)/(track->Npl()-2);
   //Defining the variable for studying the kink angle
   if (trackfill[itrk] > 1) trackfill[itrk] = 1.; //something strange happens, TO CHECK

   //Storing MC true information
   MCEventID[itrk] = track->MCEvt();
   MCTrackID[itrk] = track->MCTrack();
   MCMotherID[itrk] = track->Aid(0); //used to store MotherID information
   //number of holes and gaps            
   trk_num_holes[itrk] = track->N0();  // number of holes
   trk_max_gap[itrk] = track->CheckMaxGap();  // max of consecutive holes in a track
   
   //track angles
   TX[itrk] = track->TX();
   TY[itrk] = track->TY();
   //Information about track-vertex association
   Int_t zpos = vertexobject->GetVTa(itrk)->Zpos();
   incoming[itrk] = zpos;
   impactparameter[itrk] = vertexobject->GetVTa(itrk)->Imp();
   
   //remove tracks already associated to a vertex from the list. 
   tracklist->RemoveAt(TrackID[itrk]); //Note: since each tracks can belong to two vertex, some tracks will be tried to be removed twice, but this should be ok, RemoveAt preserves the holes
  } //end of loop on tracks
  outvertextree->Fill();   
 }//end of loop on vertices
  
 inputfile->Close();
 //Writing tree and vertex object to file
 outvertextree->Write();
 vertexrec->Write();

 savingremainingtracks(tracklist);

}

void savingremainingtracks(TObjArray* tracklist){

 //Writing list of excluded tracks
 TObjArray newtracklist = TObjArray(100000);
 int ntracks = tracklist->GetEntries();
 EdbTrackP* track = NULL;
 for (int itrk = 0; itrk < ntracks; itrk++){
  track = (EdbTrackP*) tracklist->At(itrk);
  if (track != NULL) newtracklist.Add(track);
 }

 dproc->MakeTracksTree(newtracklist, 0.,0.,"verticesandtracks.root");

}

void estimatemeanseg(EdbTrack* mytrack){ //original script by Valerio for mean seg computation
   
   EdbSegP * seg0 =0;
   EdbSegP * seg=0;
   EdbSegP * seg2=0;   
   Bool_t  * same_plate= new Bool_t[mytrack->N()]; // 
   Bool_t  * rem_seg= new Bool_t[mytrack->N()]; // 
   //Float_t  * same_plate= new Float_t[mytrack->N()]; // 
   Float_t mean_seg_x=0;
   Float_t mean_seg_y=0;
   Int_t no_same_plate=0;
   
   //Loops on segments
   
   // Find segments in the same plate
   for (int iseg = 0; iseg < mytrack->N(); iseg++){
     same_plate[iseg]=false;
     rem_seg[iseg]=false;
     seg = (EdbSegP *)(mytrack->GetSegment(iseg));
     if (iseg>0){
       seg0 = (EdbSegP *)(mytrack->GetSegment(iseg-1));
       if(seg->Plate()==seg0->Plate()){
         same_plate[iseg]=true;
         same_plate[iseg-1]=true;       
       }
     }
     //cout << "before "<< ivtx << " " << itrk << " " << iseg << " " <</* rmstransverse << " " << rmslongitudinal << " " <<*/ seg->ID()<< " " << same_plate[iseg] <<" "<<seg->X()<<" "<<seg->Y()<<" "<<seg->Z()<<" "<<seg->Plate() <<  endl;
   }
      
   // Mean x and y positions of segments in mytrack (same plate excluded)
   for (int iseg = 0; iseg < mytrack->N(); iseg++){
     seg = (EdbSegP *)(mytrack->GetSegment(iseg));
     if(same_plate[iseg]==false){
       mean_seg_x = seg->X();
       mean_seg_y = seg->Y();
       no_same_plate++;
     }
   }
   if(no_same_plate!=0){
    mean_seg_x /= no_same_plate;
    mean_seg_y /= no_same_plate;
   }
   else {
     mean_seg_x =-1;
     mean_seg_y =-1;
   }
   
   
   
   // Tag bad segments
   for (int iseg = 0; iseg < mytrack->N(); iseg++){
     seg = (EdbSegP *)(mytrack->GetSegment(iseg));
     if(iseg>0){
       seg0 = (EdbSegP *)(mytrack->GetSegment(iseg-1));
       if(seg->Plate()==seg0->Plate() && mean_seg_x!=-1 && mean_seg_y!=-1){
         float gap_seg = TMath::Sqrt(TMath::Power(seg->X()- mean_seg_x,2) + TMath::Power(seg->Y() - mean_seg_y,2)); 
         float gap_seg0 = TMath::Sqrt(TMath::Power(seg0->X()- mean_seg_x,2) + TMath::Power(seg0->Y() - mean_seg_y,2));       
         if(gap_seg > gap_seg0) rem_seg[iseg]=true;
         else rem_seg[iseg-1]=true; 
         nseg[itrk]--;                    
     }
     if(mean_seg_x==-1 && mean_seg_y==-1) {
       vertexrec->RemoveTrackFromVertex(vertexobject, itrk);
       //cout << "eccomi " << ivtx << " " << itrk << endl;
     } //rimozione di tracce inserita da Valerio
   }
   //cout << "after "<< ivtx << " " << itrk << " " << iseg << " " <</* rmstransverse << " " << rmslongitudinal << " " <<*/ seg->ID()<< " " << same_plate[iseg] << " "<< rem_seg[iseg] << " "<<seg->X()<<" "<<seg->Y()<<" "<<seg->Z()<<" "<<seg->Plate() <<  endl;
   }
}


