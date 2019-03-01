/*simple script to loop over events, MC tracks and hit points
Works without any include, provided that we are in FairShip environment,
just launch root -l simple_loop in the folder with the simulation output*/
void simple_loop(){

 TFile *file = TFile::Open("ship.conical.Pythia8CharmOnly-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> pixelpoints(reader,"PixelModulesPoint");

// cout<<"Number of events"<<reader.GetEntries()<<endl;
 int ientry = 0;
 while (reader.Next()){
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)
     cout<<tracks.GetSize()<<endl;
     //access the array of tracks
     for (const ShipMCTrack& track: tracks){
         cout<<"PdgCode of the track: "<<track.GetPdgCode()<<" at event number "<<ientry<<endl;
     }
     //access the hits:    
     for (const PixelModulesPoint& pixelpoint: pixelpoints){
        cout<<"Hit position: "<<pixelpoint.GetX()<<" at event number "<<ientry<<endl;     
     }
 } 
}