<<<<<<< Updated upstream
/*look at topology of nutau simulation for my thesis*/
void nutau_event(){
=======
/*look at topology of nutau simulation for my thesis (created by Antonio 22 June 2020)*/
Bool_t FindBrick(TGeoManager *tgeom, Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate)

void simple_loop(){
>>>>>>> Stashed changes
 //getting tree and defining arrays
 TFile *file = TFile::Open("ship.conical.Genie-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");

 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",15,0,15);
 TH2D *htauppt = new TH2D("htauppt","Transverse momentum vs momentum of tau lepton;P[GeV/c];Pt[GeV/c]",400,0,400,100,0,10);
 TH1D *htaugamma = new TH1D("htaugamma","Gamma of tau lepton;#gamma",100,0,100);
 
 int trackID, ntauhits;
 double energy, mass;

 TDatabasePDG *pdg = TDatabasePDG::Instance();

 cout<<"Number of events"<<reader.GetEntries()<<endl;
 const int nentries = reader.GetEntries();

 //vector<int> taudaughters;
 for(int ientry = 0;ientry<nentries;ientry++){
 //while (reader.Next()){
     //resetting counters
     ntauhits = 0;
     if (ientry%10000 == 0) cout<<"arrived at entry" <<ientry<<endl;
     reader.SetEntry(ientry);// keeps track of the number of event (from 0 to Nevents - 1)
     //if neutrino interacted off of our target, go to next
     if ( (TMath::Abs(tracks[0].GetStartX())>40) || (TMath::Abs(tracks[0].GetStartY())>40) ) continue;
     htauppt->Fill(tracks[1].GetP(), tracks[1].GetPt());

     //computing gamma factor E/m
     energy = tracks[1].GetEnergy();
     mass = pdg->GetParticle(tracks[1].GetPdgCode())->Mass();
     htaugamma->Fill(energy/mass);

     //access the array of tracks
/*     for (const ShipMCTrack& track: tracks){
         if (ientry < 10) cout<<"PdgCode of the track: "<<track.GetPdgCode()<<" at event number "<<ientry<<endl;
     }*/
     //access the hits: 
     
     for (const TargetPoint& targetpoint: targetpoints){
        trackID = targetpoint.GetTrackID(); 
        if (trackID == 1){ //hit from tau lepton
            ntauhits++;
        }
     } //end of hit loop
     hnfilmtau->Fill(ntauhits,tracks[0].GetWeight());
 } // end of event loop
 TCanvas *cnfilm = new TCanvas();
 hnfilmtau->Draw();

 cout<<"Fraction of short decays"<<hnfilmtau->Integral(1,1)/hnfilmtau->Integral()<<endl;
 
 TCanvas *cp = new TCanvas();
 htauppt->Draw("COLZ");

 TCanvas *cgamma = new TCanvas();
 htaugamma->Draw();
}

//find brick from position in the geometry
Bool_t FindBrick(TGeoManager *tgeom, Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate)
{
 tgeom->FindNode(x,y,z);
 if (tgeom->GetLevel() == 0) return kFALSE;//we have in the cave, no mother volume present
 const char *name = tgeom->FindNode(x,y,z)->GetMotherVolume()->GetName(); //go there
 if(strcmp(name, "Brick") == 0 ||strcmp(name, "CES") == 0){
  NPlate = tgeom->GetMother(0)->GetNumber();
  NColumn = tgeom->GetMother(2)->GetNumber();
  NRow = tgeom->GetMother(3)->GetNumber();
  NWall = tgeom->GetMother(4)->GetNumber();

  return kTRUE;
 }
 else return kFALSE;
}
