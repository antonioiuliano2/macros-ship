/*simple script to loop over events, MC tracks and hit points
Works without any include, provided that we are in FairShip environment,
just launch root -l simple_loop in the folder with the simulation output*/
#include <map>
void hit_loop(){
  std::set<int>::iterator it;

 TFile *file = TFile::Open("pythia8_Geant4_1000_0.5.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<BoxPoint> emulsionhits(reader,"BoxPoint");
 
 TH1D *hmother = new TH1D("hmother","Mother ID of tracks produced with hadronic inelastic",100,0,100);

 THStack *hsize = new THStack("hsize","Stacked multiplicity histograms");

 TH1I *hprimarysize = new TH1I("hprimarysize","Molteplicity of primary vertices",100,0,100);
 TH1I *hsecondarysize = new TH1I("hsecondarysize","Molteplicity of secondary vertices",100,0,100);
 TH1I *hvertexsize = new TH1I("hvertexsize","Molteplicity of hadron vertices",100,0,100);

 THStack *hstartz = new THStack("hstartz", "Stacked start z histograms");

 TH1D *hprimarystartz = new TH1D("hprimarystartz","Z position of primaryvertices",50,120,130);
 TH1D *hsecondarystartz = new TH1D("hsecondarystartz","Z position of secondaryvertices",50,120,130);
 TH1D *hvertexstartz = new TH1D("hvertexstartz", "Z position of vertex",50,120,130);

 TH1D *hmomentum = new TH1D("hmomentum", "momentum of tracks associated to vertex",1000,0,100);
 TH1D *hangle = new TH1D("hangle","theta angle of tracks associated to vertex", 4,-2,2);
// cout<<"Number of events"<<reader.GetEntries()<<endl;
 int ientry = 0;
 std::map<int,set<int>> vertexmap;
 while (reader.Next()){
     vertexmap.clear(); //begin of an event, clear the map
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)
     //access the array of tracks
/*     for (const ShipMCTrack& track: tracks){
         cout<<"PdgCode of the track: "<<track.GetPdgCode()<<" at event number "<<ientry<<endl;
     }*/
     //access the hits:    
     for (const BoxPoint& hit:emulsionhits){
        
        double pdgcode = hit.PdgCode();
	int trackID = hit.GetTrackID();
        double startz = tracks[trackID].GetStartZ();
        int motherID = tracks[trackID].GetMotherId();
        //Process info
        int procID = tracks[trackID].GetProcID();
        double charge;        

        if ((TDatabasePDG::Instance()->GetParticle(pdgcode))!=NULL) charge = TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge();
        else charge = 0.;

        //selection of particles to form the vertices
         
        if (procID ==23){ //Selecting hadronic inelastic interactions
        // if (p > 0.1){ //i see only tracks above 100 MeV
         hmother->Fill(motherID);         
         if (TMath::Abs(charge)>0 && tracks[trackID].GetP()>0.1) vertexmap[motherID].insert(trackID);
        // }
        }
     } //end of the loop on hits
  //readout of the vertices

  if (!vertexmap.empty()){
   //loop over the map elements
   for (std::map<int,set<int>>::iterator it=vertexmap.begin(); it!=vertexmap.end(); ++it){
    set<int> trackidset = it->second;
    if (trackidset.size()<2) continue; //at least two tracks to form a vertex
    //filling information about vertices
    hvertexsize->Fill(trackidset.size());
    hvertexstartz->Fill(tracks[(*trackidset.begin())].GetStartZ());    
  
    if (it->first == 0){ //primary proton interactions
     cout<<"primary: "<<trackidset.size()<<" "<<ientry<<endl;
     hprimarysize->Fill(trackidset.size()); //first gives the index, second the value
     hprimarystartz->Fill(tracks[(*trackidset.begin())].GetStartZ());
    }
    else{ 
     cout<<"secondary: "<<trackidset.size()<<" "<<ientry<<endl;
     hsecondarysize->Fill(trackidset.size()); //secondary hadron interactions
     hsecondarystartz->Fill(tracks[(*trackidset.begin())].GetStartZ());
    } 
   }
  }
 }//end of the loop on events 
 //***********DRAWING HISTOGRAMS****************
 hmother->Draw();

 //let's do something different, this time
 TCanvas *csize = new TCanvas();
 hsize->Add(hprimarysize);
 hsize->Add(hsecondarysize);
 hsize->Add(hvertexsize);
 hsize->Draw();
 csize->Print("molteplicity.png");

 TCanvas *cstartz = new TCanvas();
 hstartz->Add(hprimarystartz);
 hstartz->Add(hsecondarystartz);
 hstartz->Add(hvertexstartz);
 hstartz->Draw();
 cstartz->Print("startz.png");
}
