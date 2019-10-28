//read Vertex object and insert information in FairShip (created as an exercise on 23 November 2018)

#include "/afs/cern.ch/work/a/aiuliano/public/SHIPBuild_test/FairShip/shipdata/ShipMCTrack.h"
#include "/afs/cern.ch/work/a/aiuliano/public/SHIPBuild_test/FairShip/charmdet/BoxPoint.h"
#include "/cvmfs/ship.cern.ch/SHiPBuild/sw/slc6_x86-64/FairRoot/May30-ship-2/include/FairMCEventHeader.h"

void Fedra2FairShip(){

TFile *file = TFile::Open("/eos/experiment/ship/user/aiuliano/CHARM1_RUN6/Emulsion/vertexing/MVA_selection_files/vertices_firstquarter.root");
EdbVertexRec *vertexlist = (EdbVertexRec*) file->Get("EdbVertexRec"); //getting object

TFile *outfile = new TFile("/eos/experiment/ship/user/aiuliano/CHARM1_RUN6/Emulsion/vertexing/MVA_selection_files/vertices_firstquarter_formatted.root","RECREATE");
TTree *tree = new TTree("cbmsim","Reconstructed vertices, tracks and segments in emulsion");

//objects to be saved
FairMCEventHeader *savedvertex = new FairMCEventHeader(); 
TClonesArray *trackarray = new TClonesArray("ShipMCTrack");
TClonesArray *segmentarray = new TClonesArray("BoxPoint");

tree->Branch("vertex", "FairMCEventHeader",&savedvertex,32000,-1);
tree->Branch("volumetrack", &trackarray);
tree->Branch("segment", &segmentarray);
//tracks->Branch("t.","EdbSegP",&track,32000,99);
// VERTICES //

  EdbVertex *vertex = new EdbVertex();
  cout<<"number of vertices: "<<vertexlist->eVTX->GetEntries()<<endl;
  const Double_t fake_momentum = 1.0; //we do not know the real momentum, but we need a value to put the angle in MCTrack

  for(Int_t ivtx=0; ivtx<vertexlist->eVTX->GetEntries(); ivtx++){
    //clear arrays 
    trackarray->Clear("C");
    segmentarray->Clear("C");

    vertex=(EdbVertex*)(vertexlist->eVTX->At(ivtx));
    Double_t vx=vertex->X();
    Double_t vy=vertex->Y();
    Double_t vz=vertex->Z();       
    Int_t ntracks = vertex->N(); 
    cout<<"Vtx number: "<<ivtx<<" coordinates: "<<vx<<" "<<vy<<" "<<vz<<" Number of tracks: "<<ntracks<<endl;

    //Build EventHeader object
    savedvertex = new FairMCEventHeader();
    savedvertex->SetVertex(vx, vy, vz);
    savedvertex->SetNPrim(ntracks);
    cout<<"Prova Z: "<<savedvertex->GetZ()<<endl;
  //TRACKS ASSOCIATED TO VERTICES//
    for (Int_t itrk = 0; itrk < ntracks; itrk++){
     EdbTrackP *track = (EdbTrackP*) vertex->GetTrack(itrk);
     Double_t startx = track->TrackStart()->X();
     Double_t starty = track->TrackStart()->Y();
     Double_t startz = track->TrackStart()->Z();
     
     Double_t TX = track->TrackStart()->TX();
     Double_t TY = track->TrackStart()->TY();

     Double_t theta = TMath::ATan(TMath::Sqrt(TX*TX + TY*TY));
     Double_t phi = TMath::ATan2(TY,TX);
      
     Double_t px = fake_momentum * TMath::Sin(theta) * TMath::Cos(phi);
     Double_t py = fake_momentum * TMath::Sin(theta) * TMath::Sin(phi);
     Double_t pz = fake_momentum * TMath::Cos(theta);
     cout<<"Track number: "<<itrk<<" numberofsegments: "<<track->N()<<" start coordinates "<<track->TrackStart()->X()<<" "<<track->TrackStart()->Y()<<" "<<track->TrackStart()->Z()<<endl;

     //Build SHiPMCTrack object
    ShipMCTrack savedtrack  = ShipMCTrack(-1, -1, px, py,pz, 1.,startx , starty, startz, 0, track->N(), track->W()); //Constructor in the loop, DOH! Remember the ALICE boss!
    /**  Standard constructor  
    ShipMCTrack(Int_t pdgCode, Int_t motherID, Double_t px, Double_t py,
                Double_t pz, Double_t E, Double_t x, Double_t y, Double_t z,
                Double_t t, Int_t nPoints, Double_t w);**/
     cout<<"Prova startZ "<<savedtrack.GetStartZ()<<endl;
    new((*trackarray)[itrk])  ShipMCTrack( savedtrack );
    //SEGMENTS ASSOCIATED TO TRACKS
     for (Int_t iseg = 0; iseg <  track->N(); iseg++){
      EdbSegP *segment = (EdbSegP*) track->GetSegment(iseg);
      cout<<"Segment number: "<<iseg<<" position: "<<segment->X()<<" "<<segment->Y()<<" "<<segment->Z()<<endl;  

     //Build BoxPoint object
      new((*segmentarray)[iseg])  BoxPoint(itrk, segment->PID(), TVector3(segment->X(),segment->Y(),segment->Z()), TVector3(0.,0.,0.), 0.,0.,0.,0);
     } //end of segment loop
    } //end of track loop
  cout<<"End of vertex: "<<ivtx<<endl;
  cout<<endl;

  tree->Fill();
  } //end of vertex loop
 outfile->Write(); 
 vertexlist->Write();
 outfile->Close();
}
