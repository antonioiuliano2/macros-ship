/*simple script to loop over events, MC tracks and hit points
Works without any include, provided that we are in FairShip environment,
just launch root -l simple_loop in the folder with the simulation output*/

void rpc_loop(TString *inputfile);

void rpc_loop(TString *inputfile){

 TFile *file = TFile::Open(inputfile->Data()); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 //cout<<"Number of events"<<reader.GetEntries()<<endl;
 int ientry = 0;
 while (reader.Next()){
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)    
     //getting weight information
     ShipMCTrack firsttrack = tracks[0]; //weight is the same for all the tracks of the same event
     double weight = firsttrack.GetWeight(); //weight defined as sum (x_i pho_i)

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
             ShipMCTrack mytrack = tracks[trackID];
             int procID = mytrack.GetProcID();
           }
         }
     }
 } 
}

void rpc_loop(){
 TString *inputfilename = new TString("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutau/neutrinos2019/8/nu_e/ship.conical.Genie-TGeant4.root");
 rpc_loop(inputfilename);
}
