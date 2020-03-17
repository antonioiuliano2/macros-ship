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
 cout<<"Molteplicity    TrackIDmother   Ievent(startsfrom0)"<<endl;
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
        
        double px = hit.GetPx();
        double py = hit.GetPy();
        double pz = hit.GetPz();
        double momentum = pow(px*px + py*py +pz*pz,0.5);
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
         if (TMath::Abs(charge)>0 && momentum>0.1) vertexmap[motherID].insert(trackID);
        // }
        }
     } //end of the loop on hits
  //readout of the vertices

  if (!vertexmap.empty()){
   //loop over the map elements
   for (std::map<int,set<int>>::iterator it=vertexmap.begin(); it!=vertexmap.end(); ++it){
    set<int> trackidset = it->second;
    int motherpdg = tracks[it->first].GetPdgCode(); 
    double motherenergy = tracks[it->first].GetEnergy();
    if (trackidset.size()<2) continue; //at least two tracks to form a vertex
    //filling information about vertices
    hvertexsize->Fill(trackidset.size());
    hvertexstartz->Fill(tracks[(*trackidset.begin())].GetStartZ());    
  
    if ((motherpdg == 2212) && (motherenergy > 398.)){
    //if (it->first == 0){ //primary proton interactions
    // if (tracks[(*trackidset.begin())].GetStartX() < -0.25)cout<<trackidset.size()<<" "<<it->first<<" "<<ientry<<endl;
     hprimarysize->Fill(trackidset.size()); //first gives the index, second the value
     hprimarystartz->Fill(tracks[(*trackidset.begin())].GetStartZ());
    }
    else{ 
    if (tracks[(*trackidset.begin())].GetStartX() < -0.25)cout<<trackidset.size()<<" "<<it->first<<" "<<ientry<<endl;
     hsecondarysize->Fill(trackidset.size()); //secondary hadron interactions
     hsecondarystartz->Fill(tracks[(*trackidset.begin())].GetStartZ());
    } 
   }
  }
 }//end of the loop on events 
 //***********DRAWING HISTOGRAMS****************
 hmother->Draw();
 cout<<"Numero primari: "<<hprimarysize->GetEntries()<<endl;
 //let's do something different, this time
 TCanvas *csize = new TCanvas();
 hprimarysize->SetFillColor(kRed);
 hvertexsize->SetFillColor(kGreen);
 hsize->Add(hprimarysize);
 //hsize->Add(hsecondarysize);
 hsize->Add(hvertexsize);
 hsize->Draw();
 auto legend = new TLegend(0.6,0.7,0.98,0.9);
 legend->AddEntry(hprimarysize, "primary vertices");
 legend->AddEntry(hvertexsize, "all vertices, hadronic inelastic");
 legend->Draw("SAME");
 csize->Print("molteplicity.png");

 TCanvas *cstartz = new TCanvas();
 hprimarystartz->SetFillColor(kRed);
 hvertexstartz->SetFillColor(kGreen);
 hstartz->Add(hprimarystartz);
 //hstartz->Add(hsecondarystartz);
 hstartz->Add(hvertexstartz);
 hstartz->Draw();
 auto legendstartz = new TLegend(0.6,0.7,0.98,0.9);
 legendstartz->AddEntry(hprimarystartz, "primary vertices");
 legendstartz->AddEntry(hvertexstartz, "all vertices, hadronic inelastic");
 legendstartz->Draw("SAME");
 cstartz->Print("startz.png");
}
