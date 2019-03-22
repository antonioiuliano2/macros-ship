//read Vertex object and insert information in FairShip (created as an exercise on 23 November 2018)
//added MVA selection ( 22 Marzo 2019)

void save_vertices(){

TFile *file = TFile::Open("/eos/experiment/ship/user/aiuliano/CHARM1_RUN6/Emulsion/vertexing/MVA_selection_files/vertices_secondquarter.root");
EdbVertexRec *vertexlist = (EdbVertexRec*) file->Get("EdbVertexRec"); //getting object
 TTree *vtxtree = (TTree*) file->Get("vtx");

TFile *bdtfile = TFile::Open("/eos/experiment/ship/user/aiuliano/CHARM1_RUN6/Emulsion/vertexing/MVA_selection_files/vtx_BDT_evaluated.root");
TTree *bdttree = (TTree*) bdtfile->Get("bdt");

float bdt_value;
const float bdtcut = 0.; //cut from MVA discriminator

bdttree->SetBranchAddress("bdt_value",&bdt_value);

TFile *outfile = new TFile("/eos/experiment/ship/user/aiuliano/CHARM1_RUN6/Emulsion/vertexing/selected_vertices_second_quarter.root","RECREATE");
TTree *tree = new TTree("emulsion","Reconstructed vertices, tracks and segments in emulsion");

//objects to be saved
TClonesArray *trackarray = new TClonesArray("EdbSegP");
TClonesArray *segmentarray = new TClonesArray("EdbSegP");
EdbTrackP *track;
Float_t vx, vy, vz;
Int_t ntracks, vID, manualcheck;
const Int_t NTrackMax=100;
Int_t nseg[NTrackMax], trid[NTrackMax], npl[NTrackMax],n0[NTrackMax];
Int_t nsavedseg; //segments already filled

tree->Branch("vID",&vID,"vID/I"); //ID in the origina FEDRA file
tree->Branch("manualcheck",&manualcheck,"manualcheck/I");
tree->Branch("vx",&vx,"vx/F");
tree->Branch("vy",&vy,"vy/F");
tree->Branch("vz",&vz,"vz/F");
tree->Branch("ntracks",&ntracks,"ntracks/I");
// Information from good vertex tree VERTICES //
vtxtree->SetBranchAddress("vID",&vID);
vtxtree->SetBranchAddress("quality.manualcheck",&manualcheck);

//tree->Branch("volumetrack", &trackarray);
tree->Branch("volumetrack.",&trackarray);
tree->Branch("trid",&trid, "trid[ntracks]/I");
tree->Branch("nseg",&nseg,"nseg[ntracks]/I");
tree->Branch("npl",&npl,"npl[ntracks]/I");
tree->Branch("n0",&n0,"n0[ntracks]/I");

tree->Branch("segment", &segmentarray);
//tracks->Branch("t.","EdbSegP",&track,32000,99);
// VERTICES //

  EdbVertex *vertex = new EdbVertex();
//  cout<<"number of vertices: "<<vertexlist->eVTX->GetEntries()<<endl;
  cout<<"number of vertices: "<<vtxtree->GetEntries()<<endl;
  const int trmin = 4;
  const float amin = 0.01;
//  for(Int_t ivtx=0; ivtx<vertexlist->eVTX->GetEntries(); ivtx++){
  for(Int_t ivtx=0; ivtx<vtxtree->GetEntries(); ivtx++){
    //clear arrays 
    trackarray->Clear("C");
    segmentarray->Clear("C");
    nsavedseg = 0;
    vtxtree->GetEntry(ivtx);
    cout<<"PROVA:"<<vID<<endl;
    vertex=(EdbVertex*)(vertexlist->eVTX->At(vID));
//    vID = ivtx;
    vx=vertex->X();
    vy=vertex->Y();
    vz=vertex->Z();       
    ntracks = vertex->N(); 
//    cout<<"Vtx number: "<<ivtx<<" coordinates: "<<vx<<" "<<vy<<" "<<vz<<" Number of tracks: "<<ntracks<<endl;
    //*****************CUTS ON VERTEX
    //cuts on vertices    (already did on the original tree
    //if(vertex->Flag()<0)         continue;
    //if( ntracks<trmin) continue;
    //if( vertex->MaxAperture()<amin )  continue;
    bdttree->GetEntry(vID);
    if (bdt_value < bdtcut) continue;
    //*****************END OF CUTS
  //TRACKS ASSOCIATED TO VERTICES//
    for (Int_t itrk = 0; itrk < ntracks; itrk++){
     track = (EdbTrackP*) vertex->GetTrack(itrk);

     trid[itrk] = track->ID();
     nseg[itrk] = track->N();
     npl[itrk]  = track->Npl();
     n0[itrk]   = track->N0();

  //   cout<<"Track number: "<<itrk<<" numberofsegments: "<<track->N()<<" start coordinates "<<track->TrackStart()->X()<<" "<<track->TrackStart()->Y()<<" "<<track->TrackStart()->Z()<<endl;

    new((*trackarray)[itrk])  EdbSegP( *track );
    //SEGMENTS ASSOCIATED TO TRACKS
     for (Int_t iseg = 0; iseg <  nseg[itrk]; iseg++){
      EdbSegP *segment = (EdbSegP*) track->GetSegment(iseg);
    //  cout<<"Segment number: "<<iseg<<" position: "<<segment->X()<<" "<<segment->Y()<<" "<<segment->Z()<<endl;  

      new((*segmentarray)[iseg+nsavedseg])  EdbSegP( *segment); //BoxPoint(itrk, segment->PID(), TVector3(segment->X(),segment->Y(),segment->Z()), TVector3(0.,0.,0.),
     } //end of segment loop
     nsavedseg = nsavedseg + nseg[itrk];     
    } //end of track loop
//  cout<<"End of vertex: "<<ivtx<<endl;
//  cout<<endl;

  if (ntracks >= 3) tree->Fill();
  } //end of vertex loop
 outfile->Write(); 
 vertexlist->Write();
 outfile->Close();
}
