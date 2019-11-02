//read Vertex object and insert information in FairShip (created as an exercise on 23 November 2018)

#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/BoxPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/May30-ship-2/include/FairMCEventHeader.h"

void vertex_withdistributions_flags(){
 TFile *file = TFile::Open("vertices.root");
 EdbVertexRec *vertexlist = (EdbVertexRec*) file->Get("EdbVertexRec"); //getting object
 TTree *vtxtree = (TTree*)file->Get("vtx");

 TFile *outfile = new TFile("test_vertices_withflag.root","RECREATE");
 TTree *tree = new TTree("emulsion","Reconstructed vertices, tracks and segments in emulsion");


 //objects to be saved
 FairMCEventHeader *savedvertex = new FairMCEventHeader(); 
 TClonesArray *trackarray = new TClonesArray("EdbSegP");
 TClonesArray *segmentarray = new TClonesArray("EdbSegP");
 EdbTrackP *track;

 Int_t vID;
 Int_t manualcheck;

 tree->Branch("vertex", "FairMCEventHeader",&savedvertex,32000,-1);
 //tree->Branch("volumetrack", &trackarray);
 tree->Branch("segment", &segmentarray);
 tree->Branch("volumetrack","EdbSegP",&trackarray);
 tree->Branch("vID",&vID,"vID/I");
 tree->Branch("manualcheck",&manualcheck,"manualcheck/I");
 //tracks->Branch("t.","EdbSegP",&track,32000,99);
 // VERTICES //
 vtxtree->SetBranchAddress("vID",&vID);
 vtxtree->SetBranchAddress("quality.manualcheck",&manualcheck);

 EdbVertex *vertex = new EdbVertex();
 cout<<"number of vertices: "<<vtxtree->GetEntries()<<endl;


 for(Int_t ivtx=0; ivtx<vtxtree->GetEntries(); ivtx++){
    //clear arrays 
    trackarray->Clear("C");
    segmentarray->Clear("C");

    vtxtree->GetEntry(ivtx);

    vertex=(EdbVertex*)(vertexlist->eVTX->At(vID));
    Double_t vx=vertex->X();
    Double_t vy=vertex->Y();
    Double_t vz=vertex->Z();       
    Int_t ntracks = vertex->N(); 
    
    if (ntracks > 2) continue; //SELECTION OF VERTICES
    
    cout<<"Vtx number: "<<ivtx<<" coordinates: "<<vx<<" "<<vy<<" "<<vz<<" Number of tracks: "<<ntracks<<endl;

    //Build EventHeader object
    savedvertex = new FairMCEventHeader();
    savedvertex->SetVertex(vx, vy, vz);
    savedvertex->SetNPrim(ntracks);
    cout<<"Prova Z: "<<savedvertex->GetZ()<<endl;
  //TRACKS ASSOCIATED TO VERTICES//
    for (Int_t itrk = 0; itrk < ntracks; itrk++){
     track = (EdbTrackP*) vertex->GetTrack(itrk);

     cout<<"Track number: "<<itrk<<" numberofsegments: "<<track->N()<<" start coordinates "<<track->TrackStart()->X()<<" "<<track->TrackStart()->Y()<<" "<<track->TrackStart()->Z()<<endl;

     //Build SHiPMCTrack object
    //ShipMCTrack savedtrack  = ShipMCTrack(-1, -1, 0., 0.,0., 0.,startx , starty, startz, 0, track->N(), track->W()); //Constructor in the loop, DOH! Remember the ALICE boss!
    /**  Standard constructor  
    ShipMCTrack(Int_t pdgCode, Int_t motherID, Double_t px, Double_t py,
                Double_t pz, Double_t E, Double_t x, Double_t y, Double_t z,
                Double_t t, Int_t nPoints, Double_t w);**/
    // cout<<"Prova startZ "<<savedtrack.GetStartZ()<<endl;
    new((*trackarray)[itrk])  EdbSegP( *track );
    //SEGMENTS ASSOCIATED TO TRACKS
     for (Int_t iseg = 0; iseg <  track->N(); iseg++){
      EdbSegP *segment = (EdbSegP*) track->GetSegment(iseg);
      cout<<"Segment number: "<<iseg<<" position: "<<segment->X()<<" "<<segment->Y()<<" "<<segment->Z()<<endl;  

     //Build BoxPoint object
      new((*segmentarray)[iseg])  EdbSegP( *segment); //BoxPoint(itrk, segment->PID(), TVector3(segment->X(),segment->Y(),segment->Z()), TVector3(0.,0.,0.), 0.,0.,0.,0);
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

void savingvertex(){

TFile *file = TFile::Open("/home/utente/Lavoro/Analisi/CharmData/CH1-R6/vertices.root");
EdbVertexRec *vertexlist = (EdbVertexRec*) file->Get("EdbVertexRec"); //getting object

TFile *outfile = new TFile("/home/utente/Lavoro/Analisi/CharmData/CH1-R6/vertices_withtrackupstream_16squarecm.root","RECREATE");
TTree *tree = new TTree("emulsion","Reconstructed vertices, tracks and segments in emulsion");

//objects to be saved
FairMCEventHeader *savedvertex = new FairMCEventHeader(); 
TClonesArray *trackarray = new TClonesArray("EdbSegP");
TClonesArray *segmentarray = new TClonesArray("EdbSegP");
EdbTrackP *track;

tree->Branch("vertex", "FairMCEventHeader",&savedvertex,32000,-1);
//tree->Branch("volumetrack", &trackarray);
tree->Branch("segment", &segmentarray);
tree->Branch("volumetrack","EdbSegP",&trackarray);
//tracks->Branch("t.","EdbSegP",&track,32000,99);
// VERTICES //

EdbVertex *vertex = new EdbVertex();
cout<<"number of vertices: "<<vertexlist->eVTX->GetEntries()<<endl;


for(Int_t ivtx=0; ivtx<vertexlist->eVTX->GetEntries(); ivtx++){
    //clear arrays 
    trackarray->Clear("C");
    segmentarray->Clear("C");

    vertex=(EdbVertex*)(vertexlist->eVTX->At(ivtx));
    Double_t vx=vertex->X();
    Double_t vy=vertex->Y();
    Double_t vz=vertex->Z();       
    Int_t ntracks = vertex->N(); 
    
    if (ntracks > 2) continue; //SELECTION OF VERTICES
    
    cout<<"Vtx number: "<<ivtx<<" coordinates: "<<vx<<" "<<vy<<" "<<vz<<" Number of tracks: "<<ntracks<<endl;

    //Build EventHeader object
    savedvertex = new FairMCEventHeader();
    savedvertex->SetVertex(vx, vy, vz);
    savedvertex->SetNPrim(ntracks);
    cout<<"Prova Z: "<<savedvertex->GetZ()<<endl;
  //TRACKS ASSOCIATED TO VERTICES//
    for (Int_t itrk = 0; itrk < ntracks; itrk++){
     track = (EdbTrackP*) vertex->GetTrack(itrk);

     cout<<"Track number: "<<itrk<<" numberofsegments: "<<track->N()<<" start coordinates "<<track->TrackStart()->X()<<" "<<track->TrackStart()->Y()<<" "<<track->TrackStart()->Z()<<endl;

     //Build SHiPMCTrack object
    //ShipMCTrack savedtrack  = ShipMCTrack(-1, -1, 0., 0.,0., 0.,startx , starty, startz, 0, track->N(), track->W()); //Constructor in the loop, DOH! Remember the ALICE boss!
    /**  Standard constructor  
    ShipMCTrack(Int_t pdgCode, Int_t motherID, Double_t px, Double_t py,
                Double_t pz, Double_t E, Double_t x, Double_t y, Double_t z,
                Double_t t, Int_t nPoints, Double_t w);**/
    // cout<<"Prova startZ "<<savedtrack.GetStartZ()<<endl;
    new((*trackarray)[itrk])  EdbSegP( *track );
    //SEGMENTS ASSOCIATED TO TRACKS
     for (Int_t iseg = 0; iseg <  track->N(); iseg++){
      EdbSegP *segment = (EdbSegP*) track->GetSegment(iseg);
      cout<<"Segment number: "<<iseg<<" position: "<<segment->X()<<" "<<segment->Y()<<" "<<segment->Z()<<endl;  

     //Build BoxPoint object
      new((*segmentarray)[iseg])  EdbSegP( *segment); //BoxPoint(itrk, segment->PID(), TVector3(segment->X(),segment->Y(),segment->Z()), TVector3(0.,0.,0.), 0.,0.,0.,0);
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
