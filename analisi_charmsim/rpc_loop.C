/*simple script to loop over events, MC tracks and hit points
Works without any include, provided that we are in FairShip environment,
just launch root -l simple_loop in the folder with the simulation output*/
void rpc_loop(inputfile){

 TFile *file = TFile::Open(inputfile); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> rpcpoints(reader,"ShipRpcPoint");

// cout<<"Number of events"<<reader.GetEntries()<<endl;
 int ientry = 0;
 while (reader.Next()){
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)    
     //access the hits:    
     for (const ShipRpcPoint& rpcpoint: rpcpoints){
         double x = rpcpoint.GetX();
         double y = rpcpoint.GetY();
         double z = rpcpoint.GetZ();

         double px = rpcpoint.GetPx();
         double py = rpcpoint.GetPy();
         double pz = rpcpoint.GetPz();

         int detectorID = rpcpoint.GetDetectorID();
         int trackID = rpcpoint.GetDetectorID();
         bool interestingtrack = false;
         if (ientry < 10) cout<<"Hit X position: "<<rpcpoint.GetX()<<" at event number "<<ientry<<endl;     
         if (interestingtrack){
           if (trackID < 0) cout<<"cannot show track, too low kinetic energy" <<endl;
            else{
             ShipMCTrack *mytrack = tracks[trackID];
             int procID = mytrack->GetProcID();
           }
         }
     }
 } 
}
