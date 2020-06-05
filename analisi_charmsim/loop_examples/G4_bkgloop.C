//loop over stand-alone G4 pot background (27 Gennaio)
using namespace ROOT;
void G4_bkgloop(TString inputfile){
 TFile *file = TFile::Open(inputfile.Data()); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 //TFile * outputfile = new TFile("histos_primary_bkg.root","RECREATE");

 //defining histograms   
 TH1D * hstartz_primary = new TH1D("hstartz_primary","Startz position primaries;z[cm]",80,70,150);
 TH1D * hstartz_secondary = new TH1D("hstartz_secondary","Startz position secondaries;z[cm]",80,70,150);
 //TH1D * henergy = new TH1D("henergy","Distribution in energy;E[GeV]",400,0,400);
 TH1I * hnprongs_primary = new TH1I("hnprongs_primary", "Number of tracks at primary vertex",100,0,100);
 TH1I * hnprongs_secondary = new TH1I("hnprongs_secondary", "Number of tracks at secondary vertex",100,0,100);

 double endtarget = 150; //I am not interested in Goliath and MuonTagger interactions

 const int nentries = reader.GetEntries();

 int motherid, procid;
 int nprimaries;
 RVec<double> energy_primaries;
 double mindeltaE = 10.;
 double startz;

 cout<<"Starting loop over "<<nentries<<endl;

 struct myHadronicVertex {
        int molt;
        double vx,vy,vz;
 };

 map<int, myHadronicVertex> VertexfromID;
 vector<int> vertex_momIDs;

 TDatabasePDG *pdg = TDatabasePDG::Instance();
 for (int ientry = 0; ientry < nentries; ientry++){
  if (ientry % 10000==0) cout<<"Arrived at entry "<<ientry<<endl;
  reader.SetEntry(ientry);
  //starting loop over tracks
  nprimaries = 0;
  startz = 0.;
  //clearing containers (never forget, never worry)
  VertexfromID.clear();
  vertex_momIDs.clear();
  //energy_primaries.clear();
  int pdgcode;
  double primaryenergy = tracks[0].GetEnergy();
  bool interactingintarget = false;
  for (const ShipMCTrack &MCTrack:tracks){

   motherid = MCTrack.GetMotherId();
   procid = MCTrack.GetProcID();

   if (procid == 23){ // primary daughter selection
     startz = MCTrack.GetStartZ();
     pdgcode = MCTrack.GetPdgCode();
     int charge = 0;
     if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
     if(startz<endtarget&&TMath::Abs(charge)>0){ 
      nprimaries++;
     // energy_primaries.push_back(MCTrack.GetEnergy());
     vector<int>::iterator foundvertex = find (vertex_momIDs.begin(), vertex_momIDs.end(), motherid);
 
      if (foundvertex==vertex_momIDs.end()){ //build vertex struct
         myHadronicVertex *newvertex = new myHadronicVertex();
         VertexfromID[motherid] = *newvertex; 
  
         VertexfromID[motherid].molt = 0;      
         VertexfromID[motherid].vz = startz;  
  
         vertex_momIDs.push_back(motherid);
       }
       //if(motherid==7299 && ientry==0) cout<<"Eccomi "<<pdgcode<<endl;   
       VertexfromID[motherid].molt++;      
     }
    }
  } //end loop over tracks
  //filling vertex histograms
  for (int vid:vertex_momIDs){
    //if(VertexfromID[vid].molt == 1) cout<<vid<<" "<<ientry<<endl;    
    if(vid>0){
     if(VertexfromID[vid].molt>=2) hnprongs_secondary->Fill(VertexfromID[vid].molt);
     if(VertexfromID[vid].molt>=8) hstartz_secondary->Fill(VertexfromID[vid].vz);
    }
    else if (vid==0){
     if(VertexfromID[vid].molt>=2) hnprongs_primary->Fill(VertexfromID[vid].molt);
     if(VertexfromID[vid].molt>=8) hstartz_primary->Fill(VertexfromID[vid].vz);
    }
   }
 /* if (nprimaries > 0){
    if (primaryenergy - Max(energy_primaries)>mindeltaE){
     hnprimaries->Fill(nprimaries);
     hstartz->Fill(startz);
     for (const double energy: energy_primaries) henergy->Fill(energy);
    }
    //if (nprimaries < 3) cout<<"Suspect event "<<ientry<<endl;
   }*/
 }//end loop over events
 //outputfile->cd();
 //hstartz->Write();
 //henergy->Write();
 TCanvas *cnprongs = new TCanvas();
 hnprongs_secondary->SetLineColor(kRed);
 hnprongs_secondary->Draw();
 hnprongs_primary->Draw("SAMES");
 cnprongs->BuildLegend();

 TCanvas *cvz = new TCanvas();
 hstartz_secondary->SetLineColor(kRed);
 hstartz_secondary->Draw();
 hstartz_primary->Draw("SAMES");
 cvz->BuildLegend();
 //outputfile->Close();
}
