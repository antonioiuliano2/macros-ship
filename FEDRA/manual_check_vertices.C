
/*A tool for checking anomalies in the vertices distribution after looking at them at the event display. first_creation() creates a friend tree of the distribution tree, that can be set by the main script() manual_check_vertices(). Created by Antonio on 6 February 2019
 to create distribution tree launch 
 root -l 
 .L manual_check_vertices.C
 modify_distribution_tree()

Test commit 18 Aprile 2019
*/
TString path = "/afs/cern.ch/work/a/aiuliano/public/sim_fedra/CH1_testcharmpoints/b000001/"; 

TString inputfilename = path + "vertexing2/vertices_MC.root";

void first_creation(){ //to be launched the first time a new vertices tree is created (reset of indices)

 TFile *inputfile = TFile::Open(inputfilename.Data(),"UPDATE"); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;
 TTree *vertextree = (TTree*) inputfile->Get("vtx");
 
 int nentries = vertextree->GetEntries();
 Int_t vID;
 Int_t manualcheck;
 //creating a new tree
 TTree *output = new TTree("quality","Try to add manual check quality");
 output->Branch("manualcheck",&manualcheck,"manualcheck/I");
 
 cout<<"File di input con: "<<nentries<<" vertici "<<endl;
 for (int i = 0; i < nentries; i++){
  vertextree->GetEntry(i);
  manualcheck = 1000; //dummy value as a placeholder for the first time
  output->Fill();
 }
 vertextree->AddFriend(output);
 output->Write("",TObject::kOverwrite);
 vertextree->Write("",TObject::kOverwrite);
}

void manual_check_vertices(){
 //intervallo di entries del tree da indagare
 

 //opening the files and getting the trees
 //TFile *inputfile = TFile::Open("/ship/CHARM2018/CH1-R6/b000001/prova_vertexing/17_02_19/vertices.root","UPDATE");
 TFile *inputfile = TFile::Open(inputfilename.Data(),"UPDATE"); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;

 TTree *vertextree = (TTree*)inputfile->Get("vtx");
 TTree *qualitytree = (TTree*)inputfile->Get("quality");

 EdbVertexRec *gEVR = (EdbVertexRec*) inputfile->Get("EdbVertexRec");
 Int_t inputcheck;
 qualitytree->SetBranchAddress("manualcheck",&inputcheck);

 
 int nentries = vertextree->GetEntries();
 //asking to the user how many entries they want to check
 int firstentry = -1;
 int lastentry = -1;
 cout<<"Sketch of a tool for manual check of vertex display"<<endl;
 cout<<"The vertex tree contains " <<nentries<<" vertices. How many do you want to check? "<<endl;

 cout<<"From entry number: "<<endl;
 cin>>firstentry;
 cout<<"To entry number: "<<endl;
 cin>>lastentry;


 Int_t vID;
 Int_t manualcheck;
 vertextree->SetBranchAddress("vID",&vID); //I want to know the vID for printing
 //creating a new tree
 TTree *output = new TTree("quality","Try to add manual check quality");
 output->Branch("manualcheck",&manualcheck,"manualcheck/I");
 
 cout<<"File di input con: "<<nentries<<" vertici "<<endl;
 //try to insert the display
 const char *dsname="display-v";
 EdbDisplay   *ds=0;
 //vertex and track object for drawing
 EdbVertex *v = 0;
 EdbTrackP *t = 0;
 TObjArray *varr = new TObjArray();
 TObjArray *tarr = new TObjArray();
 
 for (int i = 0; i < nentries; i++){
  //clearing the arrays 
  varr->Clear();
  tarr->Clear();
  //getting the entries for both trees
  vertextree->GetEntry(i);
  qualitytree->GetEntry(i);
  if (i < firstentry || i > lastentry) manualcheck = inputcheck; //leave the previous value
  else { //set the value manually
   //adding vertex and tracks to arrays for drawing
   v = (EdbVertex *)(gEVR->eVTX->At(vID));
   varr->Add(v);
   for (int itrk = 0; itrk < v->N();itrk++){
    t = (EdbTrackP*) v->GetTrack(itrk);
    tarr->Add(t);
   }

   //display object. We create it only the first time
   ds = EdbDisplay::EdbDisplayExist(dsname);
   if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-40000., 0.);
   ds->SetVerRec(gEVR);
   ds->SetArrTr( tarr );
   printf("%d tracks to display\n", tarr->GetEntries() );
   ds->SetArrV( varr );
   printf("%d vertex to display\n", varr->GetEntries() );
   ds->SetDrawTracks(4);
   ds->SetDrawVertex(1);

   ds->Refresh(); //we refresh the display every time and we show it again
   ds->Draw();

   gPad->Modified(); //from ROOT forum about how to wait for user input within a canvas
   gPad->Update();
   gSystem->ProcessEvents();
  
   //we show some information and we ask the user about the vertex
   cout<<"Entry number "<<i<<" vertex ID "<<vID<<endl;
   cout<<"Last manual check reported (1000 being not yet setted):"<<inputcheck<<endl;
   cout<<"Please evaluate the vertex (-1 bad, 1 good, 0 unknown)"<<endl;
   cin >>manualcheck;

  }
  //Fill every time. Only the values between firstentry and lastentry should be changed
  output->Fill();
 }
 //we overwrite only the friend tree. We DO NOT want to change the vertex tree
 output->Write("",TObject::kOverwrite);
 inputfile->Close();
}


void modify_distribution_tree(){ //script to access the tree with the variables
 TFile *inputfile = TFile::Open(inputfilename.Data()); 
 if (inputfile == NULL) cout<<"ERROR: inputfile not found"<<endl;
 TTree *vertextree = (TTree*) inputfile->Get("vtx");
 //TTree *qualitytree = (TTree*) file->Get("quality");
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
 Float_t trk_max_delta_theta[maxdim];  // max delta theta between segments of a track
 Float_t trk_rms_delta_theta[maxdim];  // max delta theta between segments of a track

 TFile *outfile = new TFile("vertices_MC_modified.root","RECREATE");
 TTree *outvertextree = new TTree("vtx","Vertices");
 ////TTree *outqualitytree = qualitytree->CloneTree();
 //LIST OF VARIABLES WE WANT TO SAVE IN THE TREE
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
 //outvertextree->Branch("trk_max_delta_theta",&trk_max_delta_theta,"trk_max_delta_theta[n]/F");
 //outvertextree->Branch("trk_rms_delta_theta",&trk_rms_delta_theta,"trk_rms_delta_theta[n]/F");
 
 
 //we need some Edb objects to add new information
 EdbVertexRec *vertexrec = (EdbVertexRec*)inputfile->Get("EdbVertexRec");
 EdbVertex *vertexobject = 0;
 EdbTrackP *track = 0;
 //Vertex    *vt = 0;  da chiedere ad Antonio
 //Starting loop on vertex
 for (int ivtx = 0; ivtx < nvertices; ivtx++){
  vertextree->GetEntry(ivtx);
  vertexobject = (EdbVertex *)(vertexrec->eVTX->At(vID));
  vtx_max_aperture = vertexobject->MaxAperture();
  
  // vt->ndf() is number of degree of freedom
  // vt->chi2() is vertex chi2
  // vt->prob() is vertex probability
  // vt->rmsDistAngle() is average distance from vertex to tracks
  // vt->angle() is average angle between tracks
  
  //cout << "vtx_before " << ivtx << " " << n << " " << vertexobject->N() << endl;
  
  //Loop on tracks
  for (int itrk = 0; itrk < n; itrk++){
 
   track = vertexobject->GetTrack(itrk);
   nseg[itrk] = track->N();

   if (track->Npl() <=2) trackfill[itrk] = 1.;
   else trackfill[itrk] = (Float_t)(track->N()-2)/(track->Npl()-2);//computing the 'fill factor'
   //Defining the variable for studying the kink angle
   if (trackfill[itrk] > 1) trackfill[itrk] = 1.; //something strange happens, TO CHECK
   /*
   double rmsspace = 0.;
   double rmstransverse = 0.;
   double rmslongitudinal = 0.;
   int NKinkAngleUsed = 0;
   double kinkoutput = EdbEDAUtil::DTRMSTL1Kink(track, &rmsspace, &rmstransverse, &rmslongitudinal, &NKinkAngleUsed); 
   rmsthetatransverse[itrk] = rmstransverse;
   rmsthetalongitudinal[itrk] = rmslongitudinal;
   */
   TrackID[itrk] = track->Track(); //NOT ID(), which is resetted! The eTrack attribute of EdbSegP now allows association to original root tree
   //Storing MC true information
   MCEventID[itrk] = track->MCEvt();
   MCTrackID[itrk] = track->MCTrack();
   MCMotherID[itrk] = track->Aid(0); //used to store MotherID information
   //--------------------- Valerio workspace -----------------------------//
      
   EdbSegP * seg0 =0;
   EdbSegP * seg=0;
   EdbSegP * seg2=0;   
   Bool_t  * same_plate= new Bool_t[track->N()]; // 
   Bool_t  * rem_seg= new Bool_t[track->N()]; // 
   //Float_t  * same_plate= new Float_t[track->N()]; // 
   Float_t mean_seg_x=0;
   Float_t mean_seg_y=0;
   Int_t no_same_plate=0;
   
   //Loops on segments
   
   // Find segments in the same plate
   for (int iseg = 0; iseg < track->N(); iseg++){
     same_plate[iseg]=false;
     rem_seg[iseg]=false;
     seg = (EdbSegP *)(track->GetSegment(iseg));
     if (iseg>0){
       seg0 = (EdbSegP *)(track->GetSegment(iseg-1));
       if(seg->Plate()==seg0->Plate()){
         same_plate[iseg]=true;
         same_plate[iseg-1]=true;       
       }
     }
     //cout << "before "<< ivtx << " " << itrk << " " << iseg << " " <</* rmstransverse << " " << rmslongitudinal << " " <<*/ seg->ID()<< " " << same_plate[iseg] <<" "<<seg->X()<<" "<<seg->Y()<<" "<<seg->Z()<<" "<<seg->Plate() <<  endl;
   }
      
   // Mean x and y positions of segments in track (same plate excluded)
   for (int iseg = 0; iseg < track->N(); iseg++){
     seg = (EdbSegP *)(track->GetSegment(iseg));
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
   for (int iseg = 0; iseg < track->N(); iseg++){
     seg = (EdbSegP *)(track->GetSegment(iseg));
     if(iseg>0){
       seg0 = (EdbSegP *)(track->GetSegment(iseg-1));
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
   
   //cout << "vtx_after" << ivtx << " " << n << " " << vertexobject->N() << endl;
   
  //cout << "trk_after " << ivtx << " " << track->N() << endl;
   
   int seg_index=0;
   int tmp_id_seg=0;
   Float_t * seg_delta_theta = new Float_t[nseg[itrk]];
   
   
   for (int iseg = 0; iseg < track->N(); iseg++){
     if(rem_seg[iseg]==false){ //track->RemoveSegment(seg);
       seg = (EdbSegP *)(track->GetSegment(iseg));
       if(seg_index>0){
         seg0 = (EdbSegP *)(track->GetSegment(tmp_id_seg));
         seg_delta_theta[seg_index-1] = seg->DeltaTheta(seg0);
         //cout << "delta_theta " << ivtx << " " << itrk << " " << seg_index-1 << " " << seg_delta_theta[seg_index-1] << " " << seg->Plate() << " " << seg0->Plate() <<  endl; 
       }
       tmp_id_seg = iseg;
       seg_index++;
       nseg[itrk] = seg_index+2; // real number of segments
     }
    else {
     seg = (EdbSegP *)(track->GetSegment(iseg));
     if(rem_seg[iseg]==true)track->RemoveSegment(seg);
     
    }
   }
               
   trk_num_holes[itrk] = track->N0();  // number of holes
   trk_max_gap[itrk] = track->CheckMaxGap();  // max of consecutive holes in a track
   if(nseg[itrk]<2) {
     trk_max_delta_theta[itrk]=-1;
     trk_rms_delta_theta[itrk]=-1;
   }
   else {
     trk_max_delta_theta[itrk] = TMath::MaxElement(nseg[itrk],seg_delta_theta);
     trk_rms_delta_theta[itrk] = TMath::RMS(nseg[itrk],seg_delta_theta);
   }
   //cout << "trk "<< ivtx << " " << itrk << " " << trk_num_holes[itrk] << " " << trk_max_gap[itrk] << " " << trk_max_delta_theta[itrk] << endl;
   
   delete [] seg_delta_theta;
   delete [] same_plate;
   delete [] rem_seg;
   
   Int_t zpos = vertexobject->GetVTa(itrk)->Zpos();
   incoming[itrk] = zpos;
   TX[itrk] = track->TX();
   TY[itrk] = track->TY();
   impactparameter[itrk] = vertexobject->GetVTa(itrk)->Imp();
  }
  outvertextree->Fill();
 }
 
 /*for (int ivtx = 0; ivtx < nvertices; ivtx++){
   vertextree->GetEntry(ivtx);
   vertexobject = (EdbVertex *)(vertexrec->eVTX->At(vID));
   cout << "vtx_after" << ivtx << " " << n << " " << vertexobject->N() << endl;
 }*/
  
 inputfile->Close();
 //Writing tree and vertex object to file
 outvertextree->Write();
 vertexrec->Write();
 ////outqualitytree->Write();
}

void copyfriendtree(){ //copying a friend tree from a file to another

 //opening tree with new quality
 TFile *f = TFile::Open("/ship/CHARM2018/CH1-R6/b000001/prova_vertexing/18_02_19/modified_vertices.root");
 TTree *tree = (TTree*) f->Get("quality");
 //opening tree with distributions
 TTree *outqualitytree = (TTree*) tree->CloneTree();


 TFile *outfile = //TFile::Open("/ship/CHARM2018/CH1-R6/b000001/prova_vertexing/18_02_19/modified_vertices2.root","UPDATE");
 TFile::Open("/ship/CHARM2018/macros/Valerio/tmva_vertices.root","UPDATE");
 TTree *vertextree = (TTree*) outfile->Get("vtx");
 vertextree->AddFriend(outqualitytree);
 outqualitytree->Write("",TObject::kOverwrite);
 vertextree->Write("",TObject::kOverwrite);
 f->Close(); //closing the root file from 17_02
 outfile->Close();
}
