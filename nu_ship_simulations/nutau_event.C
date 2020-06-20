/*look at topology of nutau simulation for my thesis*/
void simple_loop(){
 //getting tree and defining arrays
 TFile *file = TFile::Open("ship.conical.Genie-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> targetpoints(reader,"TargetPoint");

 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",10,0,10);
 int trackID, ntauhits;
// cout<<"Number of events"<<reader.GetEntries()<<endl;
 int ientry = 0;
 while (reader.Next()){
     //resetting counters
     ntauhits = 0;
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)
     //if neutrino interacted off of our target, go to next
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
     hnfilmtau->Fill(ntauhits);
 } // end of event loop
 TCanvas *cnfilm = new TCanvas();
 hnfilmtau->Draw();
}